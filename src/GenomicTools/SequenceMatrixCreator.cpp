#include <SequenceMatrixCreator.hpp>

#include <string>
#include <stdexcept>       // std::invalid_argument, std::runtime_error
#include <utility>         // std::make_pair(), std::move()
#include <unordered_map>

#include <seqan/bed_io.h>  // BedFileIn, BedRecord
#include <seqan/seq_io.h>  // seqan::SeqFileIn
#include <dna_utility.hpp>
#include <Matrix2D.hpp>

SequenceMatrixCreator::SequenceMatrixCreator(const std::string& bed_file_path,
                                             const std::string& fasta_file_path,
                                             int from,
                                             int to)
    : MatrixCreator(bed_file_path,
                    from,
                    to),
      fasta_path(fasta_file_path),
      fasta_file()
{   seqan::BedRecord<seqan::Bed3> bed_line ;

    // compute number of regions
    this->open_bed_file() ;
    size_t n_row = 0 ;
    size_t n_col = to - from + 1 ;
    while(not seqan::atEnd(this->bed_file))
    {   seqan::readRecord(bed_line, this->bed_file) ;
        n_row++ ;
    }
    this->close_bed_file() ;

    // create the count matrix
    // init to 'N' because if a part of the matrix
    // cannot be filled, it wil contain stretches of
    // 'N'
    this->matrix = Matrix2D<int>(n_row, n_col, dna::char_to_int('N')) ;
}

SequenceMatrixCreator::~SequenceMatrixCreator()
{   this->close_fasta_file() ;
    // bed file closed in ~MatrixCreator()
}


Matrix2D<int> SequenceMatrixCreator::create_matrix()
{
    std::unordered_map<std::string,seqan::Dna5String> seq_map ;

    // read the fasta file and store all the sequences
    this->open_fasta_file() ;
    while(not seqan::atEnd(this->fasta_file))
    {   seqan::CharString record_id ;
        seqan::Dna5String record_seq ;
        seqan::readRecord(record_id, record_seq, this->fasta_file) ;
        std::string id = seqan::toCString(record_id) ;
        // store it
        if(seq_map.find(id) == seq_map.end())
        {   seq_map.insert(std::make_pair(std::move(id),
                                          std::move(record_seq))) ;
        }
        else
        {   char msg[4096] ;
            sprintf(msg, "Error! header %s found several times in %s",
                    id.c_str(), this->fasta_path.c_str()) ;
            throw std::runtime_error(msg) ;
        }
    }
    this->close_fasta_file() ;

    // fill the matrix
    this->open_bed_file() ;
    size_t i=0 ;
    seqan::BedRecord<seqan::Bed3> bed_line ;
    while(not seqan::atEnd(this->bed_file))
    {   seqan::readRecord(bed_line, this->bed_file) ;
        std::string region_chr = seqan::toCString(bed_line.ref) ;
        // get sequence [from, to)
        int region_mid = MatrixCreator::get_center_pos(bed_line) ;
        int region_start = std::max(0, region_mid + from) ;
        int region_end   = region_mid + to + 1 ;
        auto iter = seq_map.find(region_chr) ;
        if(iter == seq_map.end())
        {   char msg[4096] ;
            sprintf(msg, "Error! %s sequence cannot be found in %s",
                    region_chr.c_str(), this->fasta_path.c_str()) ;
            throw std::runtime_error(msg) ;
        }
        else
        {   // auto& seq_name = iter->first ;
            auto&      seq = iter->second ;
            for(int j_seq=region_start, j_mat=0;
                j_seq<region_end and j_seq<(int)seqan::length(iter->second);
                j_seq++, j_mat++)
            {   this->matrix(i,j_mat) = dna::char_to_int(seq[j_seq]) ; }
        }
        i++ ;
    }
    this->close_bed_file() ;
    return this->matrix ;
}

void SequenceMatrixCreator::open_fasta_file()
{   if(not seqan::open(this->fasta_file, this->fasta_path.c_str()))
    {   char msg[4096] ;
        sprintf(msg, "cannot open %s", this->fasta_path.c_str()) ;
         throw std::runtime_error(msg) ;
    }
}

void SequenceMatrixCreator::close_fasta_file()
{   seqan::close(this->fasta_file) ; }
