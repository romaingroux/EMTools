#include <vector>
#include <string>

#include <seqan/bed_io.h>  // BedFileIn

#include <MatrixCreator.hpp>
#include <GenomeRegion.hpp>
#include <Matrix2D.hpp>




MatrixCreator::MatrixCreator(const std::string& bed_file_path,
                             int from,
                             int to)
    : bed_path(bed_file_path),
      bed_file(),
      from(from),
      to(to),
      matrix()
{}

int MatrixCreator::get_center_pos(const seqan::BedRecord<seqan::Bed3>& bed_line)
{   int region_len = bed_line.endPos - bed_line.beginPos ;
    int region_mid = bed_line.beginPos + (region_len / 2) ;
    return region_mid ;
}

MatrixCreator::~MatrixCreator()
{   this->close_bed_file() ; }

void MatrixCreator::open_bed_file()
{   if(not seqan::open(this->bed_file, this->bed_path.c_str()))
    {   char msg[4096] ;
        sprintf(msg, "cannot open %s", this->bed_path.c_str()) ;
         throw std::runtime_error(msg) ;
    }
}

void MatrixCreator::close_bed_file()
{   seqan::close(this->bed_file) ; }
