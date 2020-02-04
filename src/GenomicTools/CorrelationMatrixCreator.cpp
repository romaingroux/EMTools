#include <string>
#include <vector>
#include <stdexcept>  // std::runtime_error

#include <seqan/bam_io.h>  // BamFileIn
#include <seqan/bed_io.h>  // BedFileIn

#include <CorrelationMatrixCreator.hpp>
#include <Matrix2D.hpp>


/* A lambda to sort GenomeRegion by ascending starting coordinate
 */
auto sortByStartPos = [](const GenomeRegion& r1, const GenomeRegion& r2) -> bool
{ return r1 < r2 ;
} ;

CorrelationMatrixCreator::CorrelationMatrixCreator(const std::string& bed_file_path,
                                                   const std::string& bam_file_path,
                                                   const std::string& bai_file_path,
                                                   int from,
                                                   int to,
                                                   int bin_size,
                                                   CorrelationMatrixCreator::methods method)
    : ReadMatrixCreator(bed_file_path,
                        bam_file_path,
                        bai_file_path,
                        from,
                        to,
                        bin_size,
                        method),
      target_list_fw(),
      target_list_rv()
{
    seqan::BedRecord<seqan::Bed3> bed_line ;

    // compute coordinates relative to each region
    this->compute_relative_bin_coord() ;
    size_t n_col = this->relative_bin_coord.size() ;

    // compute number of regions and get valid chromosomes names
    this->open_bed_file() ;
    this->open_bam_file() ;
    seqan::BamHeader header ;
    seqan::readHeader(header, bam_file) ;
    size_t n_row = 0 ;
    while(not seqan::atEnd(this->bed_file))
    {   seqan::readRecord(bed_line, this->bed_file) ;
        std::string chrom_name = seqan::toCString(bed_line.ref) ;
        // new chromosome
        if(this->chrom_map_names.find(chrom_name) ==
                this->chrom_map_names.end())
        {   int chrom_idx = -1 ;
            seqan::getIdByName(chrom_idx,
                               seqan::contigNamesCache(seqan::context(this->bam_file)),
                               chrom_name) ;
            this->chrom_map_names[chrom_name] = chrom_idx ;
        }
        n_row++ ;
    }
    this->close_bed_file() ;
    this->close_bam_file() ;

    // create the count matrix
    this->matrix = Matrix2D<int>(n_row, n_col, 0.) ;
    // create the region matrix
    this->matrix_bins =
        std::vector<std::vector<GenomeRegion>>
            (n_row,std::vector<GenomeRegion>(n_col)) ;
    this->open_bed_file() ;
    this->open_bam_file() ;
    size_t i = 0 ;
    while(not seqan::atEnd(this->bed_file))
    {   seqan::readRecord(bed_line, this->bed_file) ;
        // find the region limits
        std::string region_chr = seqan::toCString(bed_line.ref) ;
        // int region_len         = bed_line.endPos - bed_line.beginPos ;
        // int region_mid         = bed_line.beginPos + (region_len / 2) ;
        int region_mid = CorrelationMatrixCreator::get_center_pos(bed_line) ;

        // compute the absolute bins coordinates for this region
        // and create the bins in this region
        for(size_t j=0; j<n_col; j++)
        {   const auto& relative_coord = this->relative_bin_coord[j] ;
            this->matrix_bins[i][j] =
                    GenomeRegion(region_chr,
                                 this->chrom_map_names[region_chr],
                                 region_mid + relative_coord.first,
                                 region_mid + relative_coord.second) ;
        }
        i++ ;
    }
    this->close_bed_file() ;
    this->close_bam_file() ;
}

CorrelationMatrixCreator::~CorrelationMatrixCreator()
{   this->close_bam_file() ;
    // bed file is closed in ~MatrixCreator()
}

Matrix2D<int> CorrelationMatrixCreator::create_matrix()
{
    this->open_bam_file() ;
    this->open_bai_file() ;

    // read BAM header
    seqan::BamHeader bam_header ;
    seqan::readHeader(bam_header, this->bam_file) ;

    for(size_t i=0; i<this->matrix.get_nrow(); i++)
    {
        const auto& row = this->matrix_bins[i] ;
        GenomeRegion region(row.front().chromosome,
                            row.front().chromosome_idx,
                            row.front().start,
                            row.back().end) ;

        bool jump = this->jump_upstream(region, 600) ;
        if(not jump)
        {   continue ; }
        // read all relevant targets
        this->to_downstream_target(region) ;
        // update count matrix row
        this->update_count_matrix(i) ;
        // clean buffers
        this->clear_target_lists() ;
    }
    this->close_bam_file() ;
    return this->matrix ;
}

bool CorrelationMatrixCreator::jump_upstream(const GenomeRegion& region,
                                             int margin)
{   bool has_alignment = false ;
    int rID = -10 ;
    if(this->chrom_map_names.find(region.chromosome) !=
       this->chrom_map_names.end())
    {   rID = this->chrom_map_names[region.chromosome] ; }
    else
    {   char msg[4096] ;
        sprintf(msg, "Error! chromosome %s is not linked with a valid ID in BAM file",
                region.chromosome.c_str()) ;
        std::cerr << msg << std::endl ;
        return false ;
    }

    int start = std::max(0, region.start - margin) ;
    int end   = start + 1 ;
    bool jump = seqan::jumpToRegion(this->bam_file,
                                    has_alignment,
                                    rID,
                                    start,
                                    end,
                                    this->bai_file) ;
    return jump ;
}

void CorrelationMatrixCreator::to_downstream_target(const GenomeRegion& region)
{   if(this->method == CorrelationMatrixCreator::methods::READ or
       this->method == CorrelationMatrixCreator::methods::READ_ATAC)
    {   this->to_downstream_read(region) ; }
    else
    {   this->to_downstream_fragment(region) ; }
}

void CorrelationMatrixCreator::to_downstream_read(const GenomeRegion& region)
{   bool done = false ;

    seqan::BamAlignmentRecord record ;

    while(not seqan::atEnd(this->bam_file) and
          not done)
    {   // QC check and transform record
        seqan::readRecord(record, this->bam_file) ;
        if(not CorrelationMatrixCreator::is_good_read(record) or
           not this->is_valid_chromosome(record))
        {   continue ; }

        GenomeRegion target ;
        try
        {   if(this->method == CorrelationMatrixCreator::methods::READ)
            {   target = GenomeRegion::constructRead(record, this->bam_file) ; }
            else
            {   target = GenomeRegion::constructReadATAC(record, this->bam_file) ; }
        }
        catch(std::invalid_argument& e)
        {   // connect to cerr to write in SAM
            seqan::BamFileOut samFileOut(seqan::context(this->bam_file),
                                         std::cerr,
                                         seqan::Sam()) ;
            std::cerr << "std::invalid_argument caught! could not use "
                         "this record as read: " << std::endl ;
            writeRecord(samFileOut, record) ;
            std::cerr << "message was : " << e.what() << std::endl << std::endl ;
            continue ;
        }

        // upstream -> continue
        if(target < region)
        {   continue ; }
        // overlap -> store
        else if(target | region)
        {   if(not seqan::hasFlagRC(record))
            {   this->target_list_fw.push_back(target) ; }
            else
            {   this->target_list_rv.push_back(target) ; }
        }
        // downstream -> stop
        else
        {   done = true ; }
    }
}

void CorrelationMatrixCreator::to_downstream_fragment(const GenomeRegion& region)
{
    bool done = false ;

    seqan::BamAlignmentRecord record ;

    while(not seqan::atEnd(this->bam_file) and
          not done)
    {   // QC check and transform record
        seqan::readRecord(record, this->bam_file) ;
        if(not CorrelationMatrixCreator::is_good_pair(record) or
           not this->is_valid_chromosome(record))
        {   continue ; }

        GenomeRegion target ;
        try
        {   target = GenomeRegion::constructFragment(record, this->bam_file) ; }
        catch(std::invalid_argument& e)
        {   // connect to cerr to write in SAM
            seqan::BamFileOut samFileOut(seqan::context(this->bam_file),
                                         std::cerr,
                                         seqan::Sam()) ;
            std::cerr << "std::invalid_argument caught! could not use "
                         "this record as fragment: " << std::endl ;
            writeRecord(samFileOut, record) ;
            std::cerr << "message was : " << e.what() << std::endl << std::endl ;
            continue ;
        }

        // upstream -> continue
        if(target < region)
        {    continue ; }
        // overlap -> store
        else if(target | region)
        {   if(this->method == CorrelationMatrixCreator::methods::FRAGMENT_CENTER)
            {   target = GenomeRegion::constructFragmentCenter(record,
                                                               this->bam_file) ;
                if(target | region)
                {   this->target_list_fw.push_back(target) ; }
            }
            else
            {   this->target_list_fw.push_back(target) ; }
        }
        // downstream -> stop
        else if(target > region)
        {   // std::cerr << std::endl ;
            done = true ;
        }
    }
    // std::cerr << "to_downstream_fragment END" << std::endl ;
}

void CorrelationMatrixCreator::clear_target_lists()
{   this->target_list_fw.clear() ;
    this->target_list_rv.clear() ;
}

/*
void CorrelationMatrixCreator::remove_upstream_targets(const GenomeRegion& region)
{   // forward targets
    auto iter_fw = this->target_list_fw.cbegin() ;
    while(iter_fw != this->target_list_fw.end())
    {   // remove upstream reads
        if(*iter_fw < region)
        {   iter_fw = this->target_list_fw.erase(iter_fw) ; }
        // keep overlapping reads, don't stop here
        else if(*iter_fw | region)
        {   iter_fw++ ; }
        // stop at first read downstream
        else
        {   break ; }
    }
    // reverse targets
    auto iter_rv = this->target_list_rv.cbegin() ;
    while(iter_rv != this->target_list_rv.end())
    {   // remove upstream reads
        if(*iter_rv < region)
        {   iter_rv = this->target_list_rv.erase(iter_rv) ; }
        // keep overlapping reads
        else if(*iter_rv | region)
        {   iter_rv++ ; }
        // stop at first read downstream
        else
        {   break ; }
    }
}
*/

void CorrelationMatrixCreator::update_count_matrix(size_t row_index)
{
    // forward targets
    for(const auto& iter : this->target_list_fw)
    {   auto bin_start_end = CorrelationMatrixCreator::
                                 get_bin_indices(iter, this->matrix_bins[row_index]) ;
        for(int j=bin_start_end.first; j<bin_start_end.second; j++)
        {   this->matrix(row_index,j) +=
                         iter.overlap_len(this->matrix_bins[row_index][j]) ;
        }
    }
    // reverse targets
    for(const auto& iter : this->target_list_rv)
    {   auto bin_start_end = CorrelationMatrixCreator::
                                 get_bin_indices(iter, this->matrix_bins[row_index]) ;
        for(int j=bin_start_end.first; j<bin_start_end.second; j++)
        {   this->matrix(row_index,j) +=
                         iter.overlap_len(this->matrix_bins[row_index][j]) ;
        }
    }
}

/*
void CorrelationMatrixCreator::update_count_matrix_naive(size_t row_index)
{   // forward targets
    for(const auto& iter : target_list_fw)
    {   for(size_t j=0; j<this->matrix[0].size(); j++)
        {   this->matrix[row_index][j] +=
                    iter.overlap_len(this->matrix_bins[row_index][j]) ;
        }
    }
    // reverse targets
    for(const auto& iter : target_list_rv)
    {   for(size_t j=0; j<this->matrix[0].size(); j++)
        {   this->matrix[row_index][j] +=
                    iter.overlap_len(this->matrix_bins[row_index][j]) ;
        }
    }
}
*/
