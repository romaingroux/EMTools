#include <matrix_utility.hpp>

#include <vector>
#include <stdexcept>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>


Matrix2D<int> filter_rows(const std::vector<size_t>& indices,
                     const Matrix2D<int>& m)
{
    size_t nrow = m.get_nrow() ;
    size_t ncol = m.get_ncol() ;
    size_t nrow_filtered = indices.size() ;
    size_t nrow_final =  nrow - nrow_filtered ;

    // ensure that all indices are usable
    for(auto i : indices)
    {   // this index is too big
        if(i >= nrow)
        {   char msg[4096] ;
            sprintf(msg, "Error! Filter index bigger than matrix number of rows (%zu and %zu)!",
                    i, nrow) ;
            throw std::invalid_argument(msg) ;
        }
    }

    Matrix2D<int> m_final(nrow_final, ncol, 0) ;

    for(size_t i=0, i_idx=0, i_fin=0; i<nrow; i++)
    {   // skip that row
        if((i_idx<nrow_filtered) and
           (i == indices[i_idx]))
        {   i_idx++ ;
            continue ;
        }
        // copy that row
        for(size_t j=0; j<ncol; j++)
        {   m_final(i_fin, j) = m(i,j) ; }
        i_fin++ ;
    }
    return m_final ;
}

Matrix3D<double> filter_rows(const std::vector<size_t>& indices,
                        const Matrix3D<double>& m)
{
    size_t ndim1       = m.get_dim()[0] ;
    size_t ndim2       = m.get_dim()[1] ;
    size_t ndim3       = m.get_dim()[2] ;
    size_t ndim1_idx   = indices.size() ;
    size_t ndim1_final = ndim1 - ndim1_idx  ;

    // ensure that all indices are usable
    for(auto i : indices)
    {   // this index is too big
        if(i >= ndim1)
        {   char msg[4096] ;
            sprintf(msg, "Error! Filter index bigger than matrix number of rows (%zu and %zu)!",
                    i, ndim1) ;
            throw std::invalid_argument(msg) ;
        }
    }

    Matrix3D<double> m_final(ndim1_final,ndim2,ndim3,0) ;

    for(size_t i=0, i_idx=0, i_fin=0; i<ndim1; i++)
    {   // skip that row
        if((i_idx<ndim1_idx ) and
           (i == indices[i_idx]))
        {   i_idx++ ;
            continue ;
        }
        // copy that row
        for(size_t j=0; j<ndim2; j++)
        {   for(size_t k=0; k<ndim3; k++)
            {   m_final(i_fin,j,k) = m(i,j,k) ; }
        }
         i_fin++ ;
    }
    return m_final ;
}

Matrix4D<double> filter_rows(const std::vector<size_t>& indices,
                        const Matrix4D<double>& m)
{
    size_t ndim1          = m.get_dim()[0] ;
    size_t ndim2          = m.get_dim()[1] ;
    size_t ndim3          = m.get_dim()[2] ;
    size_t ndim4          = m.get_dim()[3] ;
    size_t ndim2_filtered = indices.size() ;
    size_t ndim1_final    = ndim1 - ndim2_filtered ;

    // ensure that all indices are usable
    for(auto i : indices)
    {   // this index is too big
        if(i >= ndim1)
        {   char msg[4096] ;
            sprintf(msg, "Error! Filter index bigger than matrix number of rows (%zu and %zu)!",
                    i, ndim1) ;
            throw std::invalid_argument(msg) ;
        }
    }

    Matrix4D<double> m_final(ndim1_final, ndim2, ndim3,0) ;

    for(size_t i=0, i_idx=0, i_fin=0; i<ndim1; i++)
    {   // skip that row
        if((i_idx<ndim2_filtered) and
           (i == indices[i_idx]))
        {   i_idx++ ;
            continue ;
        }
        // copy that row
        for(size_t j=0; j<ndim2; j++)
        {   for(size_t k=0; k<ndim3; k++)
            {   for(size_t l=0; l<ndim4; l++)
                {   m_final(i_fin,j,k,l) = m(i,j,k,l) ; }
            }
        }
         i_fin++ ;
    }
    return m_final ;
}
