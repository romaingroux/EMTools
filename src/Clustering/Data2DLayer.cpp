#include<Data2DLayer.hpp>

#include <stdexcept>     // std::invalid_argument
#include <utility>       // std::move()

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>


Data2DLayer::Data2DLayer()
    : DataLayer()
{}

Data2DLayer::Data2DLayer(const Matrix2D<int>& data,
                         size_t n_class,
                         size_t n_shift,
                         bool flip,
                         bool bckg_class)
    : DataLayer(n_class, n_shift, flip, bckg_class),
      data(data),
      n_row(this->data.get_nrow()),
      n_col(this->data.get_ncol())
{
    this->l_model = this->n_col - this->n_shift + 1 ;

    // models cannot be initialise here
    // as the number of categories depend
    // on the exact class
}

Data2DLayer::Data2DLayer(Matrix2D<int>&& data,
                         size_t n_class,
                         size_t n_shift,
                         bool flip,
                         bool bckg_class)
    : DataLayer(n_class, n_shift, flip, bckg_class),
      data(data),
      n_row(this->data.get_nrow()),
      n_col(this->data.get_ncol())
{
    this->l_model = this->n_col - this->n_shift + 1 ;

    // models cannot be initialise here
    // as the number of categories depend
    // on the exact class
}

Data2DLayer::Data2DLayer(const Matrix2D<int>& data,
                         const Matrix3D<double>& model,
                         bool flip,
                         bool bckg_class)
    : DataLayer(model, flip, bckg_class),
      data(data),
      n_row(this->data.get_nrow()),
      n_col(this->data.get_ncol())
{
    this->n_shift = this->n_col - this->l_model + 1 ;

    // check if model is not too long
    if(this->n_col < this->l_model)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! model is longer than data : %zu / %zu",
                this->l_model, this->n_col) ;
        throw std::invalid_argument(msg) ;
    }
}

Data2DLayer::Data2DLayer(Matrix2D<int>&& data,
                         Matrix3D<double>&& model,
                         bool flip,
                         bool bckg_class)
    : DataLayer(std::move(model), flip, bckg_class),
      data(data),
      n_row(this->data.get_nrow()),
      n_col(this->data.get_ncol())
{
    this->n_shift = this->n_col - this->l_model + 1 ;

    // check if model is not too long
    if(this->n_col < this->l_model)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! model is longer than data : %zu / %zu",
                this->l_model, this->n_col) ;
        throw std::invalid_argument(msg) ;
    }
}

Data2DLayer::~Data2DLayer()
{}

void Data2DLayer::check_loglikelihood_dim(const Matrix4D<double>& loglikelihood) const
{   if(loglikelihood.get_dim()[0] != this->n_row)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood matrix 1st dimension is not "
                "equal to data row number : %zu / %zu",
                loglikelihood.get_dim()[0], this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
    else if(loglikelihood.get_dim()[1] != this->n_class)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood matrix 2nd dimension is not "
                "equal to model class number : %zu / %zu",
                loglikelihood.get_dim()[1], this->n_class) ;
        throw std::invalid_argument(msg) ;
    }
    else if(loglikelihood.get_dim()[2] != this->n_shift)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood matrix 3rd dimension is not "
                "equal to model shift state number : %zu / %zu",
                loglikelihood.get_dim()[2], this->n_shift) ;
        throw std::invalid_argument(msg) ;
    }
    else if(loglikelihood.get_dim()[3] != this->n_flip)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood matrix 4th dimension is not "
                "equal to model flip state number : %zu / %zu",
                loglikelihood.get_dim()[3], this->n_flip) ;
        throw std::invalid_argument(msg) ;
    }
}

void Data2DLayer::check_loglikelihood_max_dim(const vector_d& loglikelihood_max) const
{   if(loglikelihood_max.size() != this->n_row)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood_max length is not "
                "equal to data row number : %zu / %zu",
                loglikelihood_max.size(), this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
}

void Data2DLayer::check_posterior_prob_dim(const Matrix4D<double>& posterior_prob) const
{   if(posterior_prob.get_dim()[0] != this->n_row)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! posterior_prob matrix 1st dimension is not "
                "equal to data row number : %zu / %zu",
                posterior_prob.get_dim()[0], this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
    else if(posterior_prob.get_dim()[1] != this->n_class)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! posterior_prob matrix 2nd dimension is not "
                "equal to model class number : %zu / %zu",
                posterior_prob.get_dim()[1], this->n_class) ;
        throw std::invalid_argument(msg) ;
    }
    else if(posterior_prob.get_dim()[2] != this->n_shift)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! posterior_prob matrix 3rd dimension is not "
                "equal to model shift state number : %zu / %zu",
                posterior_prob.get_dim()[2], this->n_shift) ;
        throw std::invalid_argument(msg) ;
    }
    else if(posterior_prob.get_dim()[3] != this->n_flip)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! posterior_prob matrix 4th dimension is not "
                "equal to model flip state number : %zu / %zu",
                posterior_prob.get_dim()[3], this->n_flip) ;
        throw std::invalid_argument(msg) ;
    }
}

