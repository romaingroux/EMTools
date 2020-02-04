#include <ModelComputer.hpp>
#include <Matrix2D.hpp>
#include <Matrix3D.hpp>

ModelComputer::ModelComputer()
    : data_layer(nullptr)
{}

ModelComputer::~ModelComputer()
{   if(this->data_layer != nullptr)
    {   delete this->data_layer ;
        this->data_layer = nullptr ;
    }
}

Matrix2D<double> ModelComputer::get_model() const
{   // the model
    Matrix3D<double> model = this->data_layer->get_model() ;
    size_t n_class = model.get_dim()[0] ;
    size_t l_model = model.get_dim()[1] ;
    size_t n_categ = model.get_dim()[2] ;
    // a nice representation of the model
    Matrix2D<double> model_nice(n_class*n_categ,
                                l_model) ;
    for(size_t i=0; i<n_class; i++)
    {   for(size_t j=0; j<n_categ; j++)
        {   size_t row = (i*n_categ) + j ;
            for(size_t k=0; k<l_model; k++)
            {   model_nice(row,k) = model(i,k,j) ; }
        }
    }

    return model_nice ;
}
