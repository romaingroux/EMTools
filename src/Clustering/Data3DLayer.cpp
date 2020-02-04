#include <Data3DLayer.hpp>

#include <Matrix3D.hpp>
#include <Matrix4D.hpp>

Data3DLayer::Data3DLayer()
    : DataLayer()
{}

Data3DLayer::Data3DLayer(const Matrix3D<double>& data,
                         size_t n_class,
                         size_t n_shift,
                         bool flip,
                         bool bckg_class)
    : DataLayer(n_class, n_shift, flip, bckg_class),
      data(data),
      n_dim1(this->data.get_dim()[0]),
      n_dim2(this->data.get_dim()[1]),
      n_dim3(this->data.get_dim()[2])
{
    this->l_model = this->n_dim2 - this->n_shift + 1 ;

    // models cannot be initialise here
    // as the number of categories depend
    // on the exact class
}

Data3DLayer::Data3DLayer(Matrix3D<double>&& data,
                         size_t n_class,
                         size_t n_shift,
                         bool flip,
                         bool bckg_class)
    : DataLayer(n_class, n_shift, flip, bckg_class),
      data(data),
      n_dim1(this->data.get_dim()[0]),
      n_dim2(this->data.get_dim()[1]),
      n_dim3(this->data.get_dim()[2])
{
    this->l_model = this->n_dim2 - this->n_shift + 1 ;

    // models cannot be initialise here
    // as the number of categories depend
    // on the exact class
}

Data3DLayer::Data3DLayer(const Matrix3D<double>& data,
                         const Matrix3D<double>& model,
                         bool flip,
                         bool bckg_class)
    : DataLayer(model, flip, bckg_class),
      data(data),
      n_dim1(this->data.get_dim()[0]),
      n_dim2(this->data.get_dim()[1]),
      n_dim3(this->data.get_dim()[2])
{
    this->n_shift = this->n_dim2 - this->l_model + 1 ;

    // check if model is not too long
    if(this->n_dim2 < this->l_model)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! model is longer than data : %zu / %zu",
                this->l_model, this->n_dim2) ;
        throw std::invalid_argument(msg) ;
    }
}

Data3DLayer::Data3DLayer(Matrix3D<double>&& data,
                         Matrix3D<double>&& model,
                         bool flip,
                         bool bckg_class)
    : DataLayer(std::move(model), flip, bckg_class),
      data(data),
      n_dim1(this->data.get_dim()[0]),
      n_dim2(this->data.get_dim()[1]),
      n_dim3(this->data.get_dim()[2])
{
    this->n_shift = this->n_dim2 - this->l_model + 1 ;

    // check if model is not too long
    if(this->n_dim2 < this->l_model)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! model is longer than data : %zu / %zu",
                this->l_model, this->n_dim2) ;
        throw std::invalid_argument(msg) ;
    }
}

Data3DLayer::~Data3DLayer()
{}

void Data3DLayer::check_loglikelihood_dim(const Matrix4D<double>& loglikelihood) const
{   if(loglikelihood.get_dim()[0] != this->n_dim1)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood matrix 1st dimension is not "
                "equal to data 1st dimension : %zu / %zu",
                loglikelihood.get_dim()[0], this->n_dim1) ;
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

void Data3DLayer::check_loglikelihood_max_dim(const vector_d& loglikelihood_max) const
{   if(loglikelihood_max.size() != this->n_dim1)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! loglikelihood_max length is not "
                "equal to data 1st dimension : %zu / %zu",
                loglikelihood_max.size(), this->n_dim1) ;
        throw std::invalid_argument(msg) ;
    }
}

void Data3DLayer::check_posterior_prob_dim(const Matrix4D<double>& posterior_prob) const
{   if(posterior_prob.get_dim()[0] != this->n_dim1)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! posterior_prob matrix 1st dimension is not "
                "equal to data 1st dimension : %zu / %zu",
                posterior_prob.get_dim()[0], this->n_dim1) ;
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
