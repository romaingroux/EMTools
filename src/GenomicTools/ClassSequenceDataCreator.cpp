#include <ClassSequenceDataCreator.hpp>

#include <string>
#include <vector>
#include <stdexcept>  // std::invalid_argument

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <dna_utility.hpp>  // char_to_int()



ClassSequenceDataCreator::ClassSequenceDataCreator(const std::string& bed_file_path,
                                                   const std::string& fasta_file_path,
                                                   const std::string& prob_file_path,
                                                   int from,
                                                   int to,
                                                   size_t class_k)
    : bed_file_path(bed_file_path),
      fasta_file_path(fasta_file_path),
      prob_file_path(prob_file_path),
      from(from),
      to(to),
      class_k(class_k-1)
{ ; }

ClassSequenceDataCreator::~ClassSequenceDataCreator()
{ ; }


Matrix3D<double> ClassSequenceDataCreator::create_matrix()
{
    // load posterior prob
    Matrix4D<double> prob ;
    prob.load(this->prob_file_path) ;

    size_t n_row   = prob.get_dim()[0] ;
    size_t n_class = prob.get_dim()[1] ;
    size_t n_shift = prob.get_dim()[2] ;
    size_t n_flip  = prob.get_dim()[3] ;
    bool     flip  = n_flip == 2 ;

    // compute class prob
    std::vector<double> prob_colsum(n_class, 0.) ;
    double tot = 0. ;
    for(size_t i=0; i<n_row; i++)
    {   for(size_t j=0; j<n_class; j++)
        {   for(size_t k=0; k<n_shift; k++)
            {   for(size_t l=0; l<n_flip; l++)
                {   double p = prob(i,j,k,l) ;
                    prob_colsum[j] += p ;
                    tot += p ;
                }
            }
        }
    }
    for(auto& p : prob_colsum)
    {   p /= tot ; }

    // class of interest (1 -> 0-based)
    if(this->class_k > n_class)
    {   char msg[4096] ;
        sprintf(msg, "Error! invalid class of interest (%zu and %zu found)",
                this->class_k, n_class) ;
        throw std::invalid_argument(msg) ;
    }

    // realigning the data using the same shift as the partition
    // will "trim" the matrix -> create a wider matrix such that
    // the realigned matrix will be of the desired wide (original
    // from and to parameters).;
    int ext   = n_shift - 1 ;
    int ext_r = ext / 2 ;
    int ext_l = ext - ext_r ;
    SequenceMatrixCreator mc(this->bed_file_path,
                             this->fasta_file_path,
                             this->from - ext_l,
                             this->to   + ext_r) ;
    Matrix2D<int> data = mc.create_matrix() ;

    // check the compatibility with the prob matrix
    if(data.get_nrow() != prob.get_dim()[0])
    {   char msg[4096] ;
        sprintf(msg, "Error! the number of references in the BED file is not equal to the"
                     "1st dimension of the probability matrix (%zu, %zu).",
                data.get_nrow(), prob.get_dim()[0]) ;
        throw(std::invalid_argument(msg)) ;
    }
    else if(data.get_nrow() < prob.get_dim()[1])
    {   char msg[4096] ;
        sprintf(msg, "Error! the number of references in the BED file is smaller than"
                     "the number of classes in the probability matrix (%zu, %zu).",
                data.get_nrow(), prob.get_dim()[1]) ;
        throw(std::invalid_argument(msg)) ;

    }
    else if(data.get_ncol() < prob.get_dim()[2])
    {   char msg[4096] ;
        sprintf(msg, "Error! the data matrix has fewer columns (bins) than the"
                     "the shifting freedom in the probability matrix (%zu, %zu).",
                data.get_ncol(), prob.get_dim()[2]) ;
        throw(std::invalid_argument(msg)) ;
    }

    // realign matrix
    size_t n_seq = data.get_nrow() ;
    size_t n_col2 = this->to - this->from + 1 ;
    // min value is 1e-100 to avoid 0 prob
    Matrix3D<double> data2(n_seq, n_col2, 4, 1e-100) ;
    // aggregate
    for(size_t i=0; i<n_seq; i++)
    {   // std::cerr << "i " << i << std::endl ;
        for(size_t s=0; s<n_shift; s++)
        {   // aggregate
            for(size_t j_fw=0, j_rv=n_col2-1; j_fw<n_col2; j_fw++, j_rv--)
            {   int base_fw = data(i,s+j_fw) ;
                // don't consider N's
                if(base_fw == dna::char_to_int('N'))
                {   continue ; }
                int base_rv = 4 - base_fw - 1 ;
                // --------- forward ---------
                if(true)
                {   data2(i,j_fw,base_fw) += prob(i,class_k,s,0) ; }

                // --------- reverse ---------
                if(flip)
                {   data2(i,j_rv,base_rv) += prob(i,class_k,s,1) ; }
            }
        }
    }
    // clean memory
    data = Matrix2D<int>() ;
    prob = Matrix4D<double>() ;

    // normalize
    for(size_t i=0; i<n_seq; i++)
    {   for(size_t j=0; j<n_col2; j++)
        {   double tot = 0. ;
            for(size_t k=0; k<4; k++)
            {   tot += data2(i,j,k) ; }
            for(size_t k=0; k<4; k++)
            {   data2(i,j,k) /= tot ; }
        }
    }

    return data2 ;
}
