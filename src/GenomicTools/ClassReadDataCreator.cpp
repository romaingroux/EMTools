#include <ClassReadDataCreator.hpp>

#include <string>
#include <vector>
#include <stdexcept>  // std::invalid_argument

#include <Matrix2D.hpp>
#include <Matrix4D.hpp>

template<class T>
std::ostream& operator << (std::ostream& stream,
                           const std::vector<T>& v)
{   for(const auto& x:v)
    {   stream << x << " " ; }
    return stream ;
}

ClassReadDataCreator::ClassReadDataCreator(const std::string& bed_file_path,
                                                             const std::string& bam_file_path,
                                                             const std::string& bai_file_path,
                                                             const std::string& prob_file_path,
                                                             int from,
                                                             int to,
                                                             int bin_size,
                                                             size_t class_k,
                                                             CorrelationMatrixCreator::methods method)
    : bed_file_path(bed_file_path),
      bam_file_path(bam_file_path),
      bai_file_path(bai_file_path),
      prob_file_path(prob_file_path),
      from(from),
      to(to),
      bin_size(bin_size),
      class_k(class_k-1), // (1-based -> 0-based)
      method(method)
{ ; }

ClassReadDataCreator::~ClassReadDataCreator()
{ ; }


Matrix2D<int> ClassReadDataCreator::create_matrix()
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

    // check class of interest
    if(this->class_k > n_class)
    {   char msg[4096] ;
        sprintf(msg, "Error! invalid class of interest (%zu and %zu found)",
                this->class_k, n_class) ;
        throw std::invalid_argument(msg) ;
    }

    // realigning the data using the same shift as the partition
    // will "trim" the matrix
    // -> create a wider matrix such that the realigned matrix will
    // be of the desired wide (original from and to parameters).
    int ext   = n_shift - 1 ;
    int ext_r = ext / 2 ;
    int ext_l = ext - ext_r ;
    int from2 = this->from - ext_l ;
    int to2   = this->to   + ext_r ;

    CorrelationMatrixCreator mc(this->bed_file_path,
                                this->bam_file_path,
                                this->bai_file_path,
                                from2,
                                to2,
                                this->bin_size,
                                this->method) ;
    Matrix2D<int> data2 = mc.create_matrix() ;
    size_t       n_col2 = data2.get_ncol() ;

    // check the compatibility with the prob matrix
    if(data2.get_nrow() != prob.get_dim()[0])
    {   char msg[4096] ;
        sprintf(msg, "Error! the number of references in the BED file is not equal to the"
                     "1st dimension of the probability matrix (%zu, %zu).",
                data2.get_nrow(), prob.get_dim()[0]) ;
        throw(std::invalid_argument(msg)) ;
    }
    else if(data2.get_nrow() < prob.get_dim()[1])
    {   char msg[4096] ;
        sprintf(msg, "Error! the number of references in the BED file is smaller than"
                     "the number of classes in the probability matrix (%zu, %zu).",
                data2.get_nrow(), prob.get_dim()[1]) ;
        throw(std::invalid_argument(msg)) ;

    }
    else if(data2.get_ncol() < prob.get_dim()[2])
    {   char msg[4096] ;
        sprintf(msg, "Error! the data matrix has fewer columns (bins) than the"
                     "the shifting freedom in the probability matrix (%zu, %zu).",
                data2.get_ncol(), prob.get_dim()[2]) ;
        throw(std::invalid_argument(msg)) ;
    }

    // realign matrix
    // size_t n_col1 = data2.get_ncol() ;
    size_t n_col3 = this->to - this->from + 1 ;
    Matrix2D<double> data3(n_row,
                           n_col3,
                           0.) ;
    for(size_t i=0; i<n_row; i++)
    {   for(size_t s=0; s<n_shift; s++)
        {   // --------------- forward ---------------
            int from_dat2_fw = s ;
            int to_dat2_fw   = from_dat2_fw + n_col3 - 1 ;
            for(int j_dat2_fw=from_dat2_fw, j_dat3_fw=0;
                    j_dat2_fw<=to_dat2_fw;
                    j_dat2_fw++, j_dat3_fw++)
            {   data3(i,j_dat3_fw) +=
                        (prob(i,class_k,s,0) *
                         data2(i,j_dat2_fw)) /
                        prob_colsum[class_k] ;
            }
            // --------------- reverse ---------------
            if(flip)
            {   int from_dat2_rev = n_col2 - 1 - s ;
                int to_dat2_rev   = from_dat2_rev - (n_col3 - 1) ;
                for(int j_dat2_rev=from_dat2_rev, j_dat3_fw=0;
                        j_dat2_rev >= to_dat2_rev;
                        j_dat2_rev--, j_dat3_fw++)
                {   data3(i,j_dat3_fw) +=
                            (prob(i,class_k,s,1) *
                             data2(i,j_dat2_rev)) /
                            prob_colsum[class_k] ;
                }
            }
        }
    }

    // clean memory
    data2 = Matrix2D<int>() ;
    prob = Matrix4D<double>() ;

    // convert to integer
    Matrix2D<int> data4(data3.get_nrow(),
                        data3.get_ncol()) ;
    for(size_t i=0; i< data4.get_nrow(); i++)
    {   for(size_t j=0; j<data4.get_ncol(); j++)
        {   data4(i,j) = (int)data3(i,j) ; }
    }

    // clean memory
    data3 = Matrix2D<double>() ;

    return data4 ;
}
