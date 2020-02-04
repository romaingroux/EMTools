
#include <utility>  // std::move()

#include <SequenceModelComputer.hpp>
#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <SequenceLayer.hpp>

SequenceModelComputer::SequenceModelComputer(Matrix2D<int>&& data,
                                             const Matrix4D<double>& post_prob,
                                             bool bckg_class,
                                             size_t n_threads)
    : ModelComputer(),
      threads(nullptr)
{
    // parameters
    size_t n_class = post_prob.get_dim()[1] ;
    size_t n_shift = post_prob.get_dim()[2] ;
    size_t n_flip  = post_prob.get_dim()[3] ;
    bool flip      = n_flip == 2 ;

    // the threads
    if(n_threads)
    {   this->threads = new ThreadPool(n_threads) ; }

    // the data and the model
    this->data_layer = new SequenceLayer(std::move(data),
                                         n_class,
                                         n_shift,
                                         flip,
                                         bckg_class) ;

    this->data_layer->update_model(post_prob,
                                   this->threads) ;
}

SequenceModelComputer::~SequenceModelComputer()
{   // threads
    if(this->threads != nullptr)
    {   this->threads->join() ;
        delete this->threads ;
        this->threads = nullptr ;
    }
    // data and model
    if(this->data_layer != nullptr)
    {   delete this->data_layer ;
        this->data_layer = nullptr ;
    }
}
