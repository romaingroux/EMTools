#include <ReadMatrixCreator.hpp>

std::pair<int, int> ReadMatrixCreator::get_bin_indices(const GenomeRegion& target,
                                                       const std::vector<GenomeRegion>& bins)
{   // the bin range and chromosome
    GenomeRegion range(bins.front().chromosome,
                        bins.front().chromosome_idx,
                        bins.front().start,
                        bins.back().end) ;
    // no overlap
    if(not (target | range))
    {   return std::make_pair(0,0) ; }
    // overlap
    else
    {   // target goes over all bins
        if(target.start <= range.start and
           target.end   >= range.end)
        {   return std::make_pair(0, bins.size()) ; }
        // partial overlap
        else
        {   int bin_start = -1 ;
            int bin_end   = -1 ;
            int bin_size = bins.front().end - bins.front().start ;

            // start
            if(target.start <= range.start)
            {   bin_start = 0 ; }
            else
            {   bin_start = (target.start - range.start) / bin_size ; }

            // end
            if(target.end >= range.end)
            {   bin_end = bins.size() ; }
            else
            {   bin_end = ((target.end - 1 - range.start) / bin_size) + 1 ; }
            return std::make_pair(bin_start, bin_end) ;
        }
    }
}

bool ReadMatrixCreator::is_good_read(const seqan::BamAlignmentRecord& record)
{
    if(seqan::hasFlagUnmapped(record) or // read unmapped flag
       seqan::hasFlagQCNoPass(record) or // not passing QC flag
       seqan::hasFlagDuplicate(record))  // PCR duplicate flag
    {   return false ; }
    return true ;
}

bool ReadMatrixCreator::is_good_pair(const seqan::BamAlignmentRecord& record)
{
    if((not seqan::hasFlagMultiple(record)) or // is paired flag
       (not seqan::hasFlagAllProper(record)))  // each read properly aligned flag
    {   return false ; }

    if((not seqan::hasFlagFirst(record)) or // read 1st in pair flag
        seqan::hasFlagLast(record))         // mate 1st in pair flag
    {   return false ; }

    // read info
    bool read_is_rev = seqan::hasFlagRC(record) ; // read is rev flag
    int  read_start  = record.beginPos ;
    // mate info
    bool mate_is_rev = seqan::hasFlagNextRC(record) ; // mate is rev flag
    int  mate_start    = record.pNext ;

    // qc
    if((not this->is_good_read(record)) or
       // --> -->
       (not read_is_rev and not mate_is_rev) or
       // <-- <--
       (read_is_rev and mate_is_rev) or
       // <-- --> 1/2
       ((read_is_rev and not mate_is_rev) and (read_start < mate_start)) or
       // <-- --> 2/2
       ((not read_is_rev and mate_is_rev) and (read_start > mate_start)))
    {   return false ; }
    return true ;
}

ReadMatrixCreator::ReadMatrixCreator(const std::string& bed_file_path,
                                     const std::string& bam_file_path,
                                     const std::string& bai_file_path,
                                     int from,
                                     int to,
                                     int bin_size,
                                     ReadMatrixCreator::methods method)
    : MatrixCreator(bed_file_path,
                    from,
                    to),
      bin_size(bin_size),
      method(method),
      relative_bin_coord(),
      bam_path(bam_file_path),
      bai_path(bai_file_path),
      bam_file(),
      bai_file(),
      chrom_map_names(),
      matrix_bins()

{   if(this->method != ReadMatrixCreator::methods::FRAGMENT and
       this->method != ReadMatrixCreator::methods::FRAGMENT_CENTER and
       this->method != ReadMatrixCreator::methods::READ and
       this->method != ReadMatrixCreator::methods::READ_ATAC)
    {   throw std::invalid_argument("Error! Unrecognized method!") ; }
}

ReadMatrixCreator::~ReadMatrixCreator()
{   // bed file is closed in ~MatrixCreator()
}

/* Initialize Histogram (table) */
/* The windows or bins are placed such that one window will be
   centered at pos 0 (odd window size), -0.5 even (window size).
   The whole range [$from,$to] will be shortened to an integer
   number of window sizes.

   Example: $from = -20, $to = 20, $ win =5;
   Windows: [-17,-13], [-12,-8], [-7,-3], [-2,2], [3,7], [8,12], [13,17]
   New range: $from = -17, $to =17
*/
void ReadMatrixCreator::compute_relative_bin_coord()
{
    int l5_p = 0 ;
    int l3_p = 0 ;

    /* begin (xb), end (xe), and center position (xe) of window near 0  */
    int xb = -this->bin_size/2; ;
    int xe = xb + this->bin_size - 1 ;
    // int xc = (xb + xe)/2 ; // unused

    if (this->from > xb)
    {   l5_p = (this->from - xb)/this->bin_size + 1 ; }
    else
    {   l5_p = -(xb - this->from)/this->bin_size ; }
    if (this->to >= xe)
    {   l3_p = (this->to - xe)/this->bin_size ; }
    else
    {   l3_p = -(xe - this->to)/this->bin_size + 1 ; }

    /* New range  */
    this->from = xb + l5_p * this->bin_size;
    this->to = xe + l3_p * this->bin_size;

    // contains the bin coordinate limits [from,to)
    // from is the 1st position within the bin and to the
    // first position after the bin.
    size_t n_bin = ((this->to-this->from)/this->bin_size) + 1 ;

    this->relative_bin_coord = std::vector<std::pair<int,int>>(n_bin) ;

    int inf = this->from ;
    int sup = inf + this->bin_size - 1 ;
    for(size_t i=0; inf<=to; inf+=this->bin_size, sup+=this->bin_size, i++)
    {   this->relative_bin_coord[i] = std::make_pair(inf, sup+1) ; }
}

bool ReadMatrixCreator::is_valid_chromosome(const seqan::BamAlignmentRecord& record)
{
    std::string name = seqan::toCString(
                          seqan::getContigName(
                            record, this->bam_file)) ;

    if(this->chrom_map_names.find(name) == this->chrom_map_names.end())
    {   return false ; }
    return true ;
}

void ReadMatrixCreator::open_bam_file()
{    if(not seqan::open(this->bam_file, this->bam_path.c_str()))
     {   char msg[4096] ;
         sprintf(msg, "cannot open %s", this->bam_path.c_str()) ;
         throw std::runtime_error(msg) ;
     }
}

void ReadMatrixCreator::open_bai_file()
{   if(not seqan::open(this->bai_file, this->bai_path.c_str()))
    {   char msg[4096] ;
        sprintf(msg, "cannot open %s", this->bai_path.c_str()) ;
        throw std::runtime_error(msg) ;
    }
}

void ReadMatrixCreator::close_bam_file()
{   seqan::close(this->bam_file) ; }
