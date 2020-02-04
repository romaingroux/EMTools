#include <DataLayer.hpp>

#include <stdexcept>  // std::invalid_argument
#include <cmath>      // log()

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

DataLayer::DataLayer()
{ ; }

DataLayer::DataLayer(size_t n_class,
                     size_t n_shift,
                     bool flip,
                     bool bckg_class)
    : flip(flip),
      n_class(n_class),
      n_shift(n_shift),
      n_flip(flip + 1),
      bckg_class(bckg_class)
{   // models cannot be initialise here
    // as the number of categories depend
    // on the exact class
}

DataLayer::DataLayer(const Matrix3D<double>& model,
                     bool flip,
                     bool bckg_class)
    : model(model),
      flip(flip),
      n_class(this->model.get_dim()[0]),
      l_model(this->model.get_dim()[1]),
      n_category(this->model.get_dim()[2]),
      n_flip(flip + 1),
      bckg_class(bckg_class)
{ ; }

DataLayer::DataLayer(Matrix3D<double>&& model,
                     bool flip,
                     bool bckg_class)
    : model(model),
      flip(flip),
      n_class(this->model.get_dim()[0]),
      l_model(this->model.get_dim()[1]),
      n_category(this->model.get_dim()[2]),
      n_flip(flip + 1),
      bckg_class(bckg_class)
{ ; }

DataLayer::~DataLayer()
{ ; }

Matrix3D<double> DataLayer::get_model() const
{   return this->model ; }


const double DataLayer::p_min = 1e-300 ;
const double DataLayer::p_min_log = log(DataLayer::p_min) ;
