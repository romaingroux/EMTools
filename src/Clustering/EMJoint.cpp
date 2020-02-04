
#include <EMJoint.hpp>

#include <string>
#include <vector>
#include <stdexcept>
#include <future>        // std::promise, std::future
#include <utility>       // std::pair, std::move()
#include <functional>    // std::bind(), std::ref()
#include <cmath>         // exp()

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ReadLayer.hpp>
#include <SequenceLayer.hpp>
#include <ThreadPool.hpp>
#include <RandomNumberGenerator.hpp> // getRandomNumberGenerator()
#include <ConsoleProgressBar.hpp>    // ConsoleProgressBar

template<class T>
std::ostream& operator << (std::ostream& stream, const std::vector<T>& v)
{   for(const auto& x : v)
    {   stream << x << " " ; }
    return stream ;
}

EMJoint::EMJoint(const std::vector<Matrix2D<int>>& read_matrices,
                 size_t n_class,
                 size_t n_iter,
                 size_t n_shift,
                 bool flip,
                 bool bckg_class,
                 const std::string& seed,
                 size_t n_threads)
    : EMBase(read_matrices[0].get_nrow(),
             read_matrices[0].get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      n_layer(read_matrices.size()),
      loglikelihood_layer(n_layer,
                          Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.)),
      loglikelihood_max(this->n_layer,
                        vector_d(this->n_row, 0.)),
      read_layers(),
      seq_layer(nullptr),
      consseq_layer(nullptr)
{
    // check data matrices and their dimensions
    if(this->n_layer == 0)
    {   throw std::invalid_argument("Error! No data layer given!") ; }
    for(const auto& matrix : read_matrices)
    {   if(matrix.get_nrow() != this->n_row)
        {   char msg[4096] ;
            sprintf(msg, "Error! Read layers have variable row numbers "
                         "(%zu and %zu)!",
                    matrix.get_nrow(), this->n_row) ;
            throw std::invalid_argument(msg) ;
        }
        else if(matrix.get_ncol() != this->n_col)
        {   char msg[4096] ;
            sprintf(msg, "Error! Read layers have variable column numbers "
                         "(%zu and %zu)!",
                    matrix.get_ncol(), this->n_col) ;
            throw std::invalid_argument(msg) ;
        }
    }

    // initialise post prob randomly
    this->set_post_prob_random(seed) ;

    // data and models
    // create read layer and initialise the models from the post prob
    for(auto& matrix : read_matrices)
    {   // create the layer
        this->read_layers.push_back(new ReadLayer(matrix,
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class,
                                                  this->threads)) ;
        this->read_layers.back()->update_model(this->post_prob,
                                               this->threads) ;
    }
}

EMJoint::EMJoint(std::vector<Matrix2D<int>>&& read_matrices,
                 size_t n_class,
                 size_t n_iter,
                 size_t n_shift,
                 bool flip,
                 bool bckg_class,
                 const std::string& seed,
                 size_t n_threads)
    : EMBase(read_matrices[0].get_nrow(),
             read_matrices[0].get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      n_layer(read_matrices.size()),
      loglikelihood_layer(n_layer,
                          Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.)),
      loglikelihood_max(this->n_layer,
                        vector_d(this->n_row, 0.)),
      read_layers(),
      seq_layer(nullptr),
      consseq_layer(nullptr)
{
    // check data matrices and their dimensions
    if(this->n_layer == 0)
    {   throw std::invalid_argument("Error! No data layer given!") ; }
    for(const auto& matrix : read_matrices)
    {   if(matrix.get_nrow() != this->n_row)
        {   char msg[4096] ;
            sprintf(msg, "Error! Read layers have variable row numbers "
                         "(%zu and %zu)!",
                    matrix.get_nrow(), this->n_row) ;
            throw std::invalid_argument(msg) ;
        }
        else if(matrix.get_ncol() != this->n_col)
        {   char msg[4096] ;
            sprintf(msg, "Error! Read layers have variable column numbers "
                         "(%zu and %zu)!",
                    matrix.get_ncol(), this->n_col) ;
            throw std::invalid_argument(msg) ;
        }
    }

    // initialise post prob randomly
    this->set_post_prob_random(seed) ;

    // data and models
    // create read layer and initialise the models from the post prob
    for(auto& matrix : read_matrices)
    {
        // create the layer
        this->read_layers.push_back(new ReadLayer(std::move(matrix),
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class,
                                                  this->threads)) ;
        this->read_layers.back()->update_model(this->post_prob,
                                               this->threads) ;
    }
}

EMJoint::EMJoint(const std::vector<Matrix2D<int>>& read_matrices,
                 const Matrix2D<int>& seq_matrix,
                 size_t n_class,
                 size_t n_iter,
                 size_t n_shift,
                 bool flip,
                 bool bckg_class,
                 const std::string& seed,
                 size_t n_threads)
    : EMBase(read_matrices[0].get_nrow(),
             read_matrices[0].get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      n_layer(read_matrices.size()+1),
      loglikelihood_layer(this->n_layer,
                          Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.)),
      loglikelihood_max(this->n_layer,
                        vector_d(this->n_row, 0.)),
      read_layers(),
      seq_layer(nullptr),
      consseq_layer(nullptr)
{   // check data matrices and their dimensions
    for(const auto& matrix : read_matrices)
    {   if(matrix.get_nrow() != this->n_row)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix row number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_nrow(), this->n_row) ;
            throw std::invalid_argument(msg) ;
        }
        else if(matrix.get_ncol() != this->n_col)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix column number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_ncol(), this->n_col) ;
            throw std::invalid_argument(msg) ;
        }
    }
    if(seq_matrix.get_nrow() != this->n_row)
    {   char msg[4096] ;
        sprintf(msg, "Error! The sequence matrix row number is different than expected "
                     "(%zu instead of %zu)!",
                seq_matrix.get_nrow(), this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
    else if(seq_matrix.get_ncol() != this->n_col)
    {   char msg[4096] ;
        sprintf(msg, "Error! The sequence matrix column number is different than expected "
                     "(%zu instead of %zu)!",
                seq_matrix.get_ncol(), this->n_col) ;
        throw std::invalid_argument(msg) ;
    }

    // initialise post prob randomly
    // getRandomGenerator(seed) ;
    this->set_post_prob_random(seed) ;

    // data and models
    // create read layer and initialise the models from the post prob
    for(auto& matrix : read_matrices)
    {   // create the layer
        this->read_layers.push_back(new ReadLayer(matrix,
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class,
                                                  this->threads)) ;
        this->read_layers.back()->update_model(this->post_prob,
                                               this->threads) ;
    }
    // create sequence layer and initialise the models from the post prob
    this->seq_layer = new SequenceLayer(seq_matrix,
                                        this->n_class,
                                        this->n_shift,
                                        this->flip,
                                        bckg_class) ;
    this->seq_layer->update_model(this->post_prob,
                                  this->threads) ;
}

EMJoint::EMJoint(std::vector<Matrix2D<int>>&& read_matrices,
                 Matrix2D<int>&& seq_matrix,
                 size_t n_class,
                 size_t n_iter,
                 size_t n_shift,
                 bool flip,
                 bool bckg_class,
                 const std::string& seed,
                 size_t n_threads)
    : EMBase(read_matrices[0].get_nrow(),
             read_matrices[0].get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      n_layer(read_matrices.size()+1),
      loglikelihood_layer(this->n_layer,
                          Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.)),
      loglikelihood_max(this->n_layer,
                        vector_d(this->n_row, 0.)),
      read_layers(),
      seq_layer(nullptr),
      consseq_layer(nullptr)
{   // check data matrices and their dimensions
    for(const auto& matrix : read_matrices)
    {   if(matrix.get_nrow() != this->n_row)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix row number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_nrow(), this->n_row) ;
            throw std::invalid_argument(msg) ;
        }
        else if(matrix.get_ncol() != this->n_col)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix column number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_ncol(), this->n_col) ;
            throw std::invalid_argument(msg) ;
        }
    }
    if(seq_matrix.get_nrow() != this->n_row)
    {   char msg[4096] ;
        sprintf(msg, "Error! The sequence matrix row number is different than expected "
                     "(%zu instead of %zu)!",
                seq_matrix.get_nrow(), this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
    else if(seq_matrix.get_ncol() != this->n_col)
    {   char msg[4096] ;
        sprintf(msg, "Error! The sequence matrix column number is different than expected "
                     "(%zu instead of %zu)!",
                seq_matrix.get_ncol(), this->n_col) ;
        throw std::invalid_argument(msg) ;
    }

    // initialise post prob randomly
    // getRandomGenerator(seed) ;
    this->set_post_prob_random(seed) ;

    // data and models
    // create read layer and initialise the models from the post prob
    for(auto& matrix : read_matrices)
    {
        // create the layer
        this->read_layers.push_back(new ReadLayer(std::move(matrix),
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class,
                                                  this->threads)) ;
        this->read_layers.back()->update_model(this->post_prob,
                                               this->threads) ;
    }

    // create sequence layer and initialise the models from the post prob
    this->seq_layer = new SequenceLayer(std::move(seq_matrix),
                                        this->n_class,
                                        this->n_shift,
                                        this->flip,
                                        bckg_class) ;
    // intialise the models with the post prob
    this->seq_layer->update_model(this->post_prob,
                                  this->threads) ;
}

EMJoint::EMJoint(const std::vector<Matrix2D<int>>& read_matrices,
                 const Matrix3D<double>& consseq_matrix,
                 size_t n_class,
                 size_t n_iter,
                 size_t n_shift,
                 bool flip,
                 bool bckg_class,
                 const std::string& seed,
                 size_t n_threads)
    : EMBase(read_matrices[0].get_nrow(),
             read_matrices[0].get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      n_layer(read_matrices.size()+1),
      loglikelihood_layer(this->n_layer,
                          Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.)),
      loglikelihood_max(this->n_layer,
                        vector_d(this->n_row, 0.)),
      read_layers(),
      seq_layer(nullptr),
      consseq_layer(nullptr)
{   // check data matrices and their dimensions
    for(const auto& matrix : read_matrices)
    {   if(matrix.get_nrow() != this->n_row)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix row number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_nrow(), this->n_row) ;
            throw std::invalid_argument(msg) ;
        }
        else if(matrix.get_ncol() != this->n_col)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix column number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_ncol(), this->n_col) ;
            throw std::invalid_argument(msg) ;
        }
    }
    if(consseq_matrix.get_dim()[0] != this->n_row)
    {   char msg[4096] ;
        sprintf(msg, "Error! The consensus sequence matrix row number is different than expected "
                     "(%zu instead of %zu)!",
                consseq_matrix.get_dim()[0], this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
    else if(consseq_matrix.get_dim()[1] != this->n_col)
    {   char msg[4096] ;
        sprintf(msg, "Error! The sequence matrix column number is different than expected "
                     "(%zu instead of %zu)!",
                consseq_matrix.get_dim()[0], this->n_col) ;
        throw std::invalid_argument(msg) ;
    }

    // initialise post prob randomly
    // getRandomGenerator(seed) ;
    this->set_post_prob_random(seed) ;

    // data and models
    // create read layer and initialise the models from the post prob
    for(auto& matrix : read_matrices)
    {   // create the layer
        this->read_layers.push_back(new ReadLayer(matrix,
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class,
                                                  this->threads)) ;
        this->read_layers.back()->update_model(this->post_prob,
                                               this->threads) ;
    }
    // create consensus sequence layer and initialise the models from the post prob
    {
        this->consseq_layer = new ConsensusSequenceLayer(consseq_matrix,
                                                         this->n_class,
                                                         this->n_shift,
                                                         this->flip,
                                                         bckg_class) ;
        this->consseq_layer->update_model(this->post_prob,
                                          this->threads) ;
    }
}

EMJoint::EMJoint(std::vector<Matrix2D<int>>&& read_matrices,
                 Matrix3D<double>&& consseq_matrix,
                 size_t n_class,
                 size_t n_iter,
                 size_t n_shift,
                 bool flip,
                 bool bckg_class,
                 const std::string& seed,
                 size_t n_threads)
    : EMBase(read_matrices[0].get_nrow(),
             read_matrices[0].get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      n_layer(read_matrices.size()+1),
      loglikelihood_layer(this->n_layer,
                          Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.)),
      loglikelihood_max(this->n_layer,
                        vector_d(this->n_row, 0.)),
      read_layers(),
      seq_layer(nullptr),
      consseq_layer(nullptr)
{   // check data matrices and their dimensions
    for(const auto& matrix : read_matrices)
    {   if(matrix.get_nrow() != this->n_row)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix row number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_nrow(), this->n_row) ;
            throw std::invalid_argument(msg) ;
        }
        else if(matrix.get_ncol() != this->n_col)
        {   char msg[4096] ;
            sprintf(msg, "Error! A read matrix column number is different than expected "
                         "(%zu instead of %zu)!",
                    matrix.get_ncol(), this->n_col) ;
            throw std::invalid_argument(msg) ;
        }
    }
    if(consseq_matrix.get_dim()[0] != this->n_row)
    {   char msg[4096] ;
        sprintf(msg, "Error! The consensus sequence matrix row number is different than expected "
                     "(%zu instead of %zu)!",
                consseq_matrix.get_dim()[0], this->n_row) ;
        throw std::invalid_argument(msg) ;
    }
    else if(consseq_matrix.get_dim()[1] != this->n_col)
    {   char msg[4096] ;
        sprintf(msg, "Error! The sequence matrix column number is different than expected "
                     "(%zu instead of %zu)!",
                consseq_matrix.get_dim()[0], this->n_col) ;
        throw std::invalid_argument(msg) ;
    }

    // initialise post prob randomly
    // getRandomGenerator(seed) ;
    this->set_post_prob_random(seed) ;

    // data and models
    // create read layer and initialise the models from the post prob
    for(auto& matrix : read_matrices)
    {   // create the layer
        this->read_layers.push_back(new ReadLayer(std::move(matrix),
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class,
                                                  this->threads)) ;
        this->read_layers.back()->update_model(this->post_prob,
                                               this->threads) ;
    }
    // create consensus sequence layer and initialise the models from the post prob
    {
        this->consseq_layer = new ConsensusSequenceLayer(std::move(consseq_matrix),
                                                         this->n_class,
                                                         this->n_shift,
                                                         this->flip,
                                                         bckg_class) ;
        this->consseq_layer->update_model(this->post_prob,
                                          this->threads) ;
    }
}

EMJoint::~EMJoint()
{   // join the threads in case
    // deleted by EMBase destructor
    if(this->threads != nullptr)
    {   this->threads->join() ;
        delete this->threads ;
        this->threads = nullptr ;
    }

    // read data and models
    for(auto& ptr : this->read_layers)
    {   if(ptr != nullptr)
        {   delete ptr ;
            ptr = nullptr ;
        }
    }
    // sequence data and models
    if(this->seq_layer != nullptr)
    {   delete seq_layer ;
        this->seq_layer = nullptr ;
    }
    // consensus sequence data and models
    if(this->consseq_layer != nullptr)
    {   delete this->consseq_layer ;
        this->consseq_layer = nullptr ;
    }
}

std::vector<Matrix3D<double>> EMJoint::get_read_models() const
{   std::vector<Matrix3D<double>> models ;
    for(const auto& ptr : this->read_layers)
    {   models.push_back(ptr->get_model()) ; }
    return models ;
}

Matrix3D<double> EMJoint::get_sequence_models() const
{   if(this->seq_layer != nullptr)
    {   return this->seq_layer->get_model() ; }
    return Matrix3D<double>() ;
}

Matrix3D<double> EMJoint::get_consensus_sequence_models() const
{   if(this->consseq_layer != nullptr)
    {   return this->consseq_layer->get_model() ; }
    return Matrix3D<double>() ;
}


EMJoint::exit_codes EMJoint::classify()
{   size_t bar_update_n = this->n_iter ;
    ConsoleProgressBar bar(std::cerr, bar_update_n, 60, "classifying") ;


    // optimize the partition
    for(size_t n_iter=0; n_iter<this->n_iter; n_iter++)
    {   // E-step
        this->compute_loglikelihood() ;
        this->compute_post_prob() ;
        // M-step
        this->compute_class_prob() ;
        this->update_models() ;
        this->center_post_state_prob() ;
        bar.update() ;
    }
    bar.update() ; std::cerr << std::endl ;

    return EMJoint::exit_codes::ITER_MAX ;
}

void EMJoint::compute_loglikelihood()
{   // compute the loglikelihood for each layer
    size_t i = 0 ;
    for(auto& ptr : this->read_layers)
    {   ptr->compute_loglikelihoods(this->loglikelihood_layer[i],
                                    this->loglikelihood_max[i],
                                    this->threads) ;
        i++ ;
    }
    if(this->seq_layer != nullptr)
    {   this->seq_layer->compute_loglikelihoods(this->loglikelihood_layer[i],
                                                this->loglikelihood_max[i],
                                                this->threads) ;
        i++ ;
    }
    if(this->consseq_layer != nullptr)
    {   this->consseq_layer->compute_loglikelihoods(this->loglikelihood_layer[i],
                                                    this->loglikelihood_max[i],
                                                    this->threads) ;
        i++ ;
    }

    // sum the likelihood for each state, over all layers
    // and rescale the values
    // don't parallelize
    if(this->threads == nullptr)
    {   std::promise<bool> promise ;
        std::future<bool> future = promise.get_future() ;
        this->compute_loglikelihood_routine(0,
                                            this->n_row,
                                            promise) ;
        future.get() ;
    }
    // parallelize
    else
    {   size_t n_threads = this->threads->getNThread() ;

        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_row,n_threads) ;

        // get promises and futures
        std::vector<std::promise<bool>> promises(n_threads) ;
        std::vector<std::future<bool>>  futures(n_threads) ;
        for(size_t j=0; j<n_threads; j++)
        {   futures[j] = promises[j].get_future() ; }

        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t j=0; j<n_threads; j++)
        {   auto slice = slices[j] ;
            this->threads->addJob(std::move(
                                      std::bind(&EMJoint::compute_loglikelihood_routine,
                                                this,
                                                slice.first,
                                                slice.second,
                                                std::ref(promises[j])))) ;
        }
        // wait until all threads are done working
        for(auto& future : futures)
        {   future.get() ; }
        // -------------------------- threads stop ---------------------------
    }
}

void EMJoint::compute_loglikelihood_routine(size_t from,
                                            size_t to,
                                            std::promise<bool>& done)
{   // the max likelihood found per row
    std::vector<double> rowmax(to-from, std::numeric_limits<double>::lowest()) ;

    // sum over layers
    for(size_t i=from, i_rowmax=0; i<to; i++, i_rowmax++)
    {   for(size_t j=0; j<this->n_class; j++)
        {   for(size_t k=0; k<this->n_shift; k++)
            {   for(size_t l=0; l<this->n_flip; l++)
                {
                    // reset
                    this->loglikelihood(i,j,k,l) = 0. ;
                    // sum
                    for(size_t m=0; m<this->n_layer; m++)
                    {   // add rescaled layer value
                        this->loglikelihood(i,j,k,l) +=
                                (this->loglikelihood_layer[m](i,j,k,l) -
                                 this->loglikelihood_max[m][i]) ;
                    }
                    // keep track of max
                    if(this->loglikelihood(i,j,k,l) > rowmax[i_rowmax])
                    {   rowmax[i_rowmax] = this->loglikelihood(i,j,k,l) ; }
                }
            }
        }
    }

    // rescale
    for(size_t i=from, i_rowmax=0; i<to; i++, i_rowmax++)
    {   for(size_t j=0; j<this->n_class; j++)
        {   for(size_t k=0; k<this->n_shift; k++)
            {   for(size_t l=0; l<this->n_flip; l++)
                {   // this->loglikelihood(i,j,k,l) -= rowmax[i_rowmax] ;
                    this->loglikelihood(i,j,k,l) =
                            std::max(this->loglikelihood(i,j,k,l) - rowmax[i_rowmax],
                                     ReadLayer::p_min_log) ;
                }
            }
        }
    }
    done.set_value(true) ;
}


void EMJoint::compute_post_prob()
{   // don't parallelize
    if(this->threads == nullptr)
    {   std::promise<vector_d> promise ;
        std::future<vector_d> future = promise.get_future() ;
        this->compute_post_prob_routine(0, this->n_row, promise) ;
        // compute the sum of post prob and the per class sum of post prob
        // from the partial results computed on each slice
        this->post_prob_tot = 0. ;
        this->post_prob_colsum = future.get() ;
        for(const auto& prob : this->post_prob_colsum)
        {   this->post_prob_tot += prob ; }
    }
    // parallelize
    else
    {    size_t n_threads = this->threads->getNThread() ;

        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_row,n_threads) ;

        // get promises and futures
        // the function run by the threads will compute
        // the partial sum per class of post_prob for the given slice
        // this should be used to compute the complete sum of post_prob
        // and the complete sum per class of post_prob
        std::vector<std::promise<vector_d>> promises(n_threads) ;
        std::vector<std::future<vector_d>>  futures(n_threads) ;
        for(size_t i=0; i<n_threads; i++)
        {   futures[i] = promises[i].get_future() ; }

        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t i=0; i<n_threads; i++)
        {   auto slice = slices[i] ;
            this->threads->addJob(std::move(
                                      std::bind(&EMJoint::compute_post_prob_routine,
                                                this,
                                                slice.first,
                                                slice.second,
                                                std::ref(promises[i])))) ;
        }
        // wait until all threads are done working
        // compute the sum of post prob and the per class sum of post prob
        // from the partial results computed on each slice
        this->post_prob_tot = 0. ;
        this->post_prob_colsum = vector_d(this->n_class, 0.) ;
        for(auto& future : futures)
        {   auto probs = future.get() ;
            for(size_t i=0; i<this->n_class; i++)
            {   double prob = probs[i] ;
                this->post_prob_colsum[i] += prob ;
                this->post_prob_tot       += prob ;
            }
        }
        // -------------------------- threads stop ---------------------------
    }
}

void EMJoint::compute_post_prob_routine(size_t from,
                                        size_t to,
                                        std::promise<vector_d>& post_prob_colsum)
{   vector_d colsums(this->n_class, 0.) ;

    // post prob
    for(size_t i=from; i<to; i++)
    {   // reset row sum to 0
        this->post_prob_rowsum[i] = 0. ;
        for(size_t n_class=0; n_class<this->n_class; n_class++)
        {   for(size_t n_shift=0; n_shift<this->n_shift; n_shift++)
            {   for(size_t n_flip=0; n_flip<this->n_flip; n_flip++)
                {
                    double p = exp(this->loglikelihood(i,n_class,n_shift,n_flip)) *
                              this->post_state_prob(n_class,n_shift,n_flip) ;
                    // double p = std::max(exp(this->loglikelihood(i,n_class,n_shift,n_flip)) *
                    //                     this->post_state_prob(n_class,n_shift,n_flip),
                    //                     ReadLayer::p_min) ;
                    this->post_prob(i,n_class,n_shift,n_flip) = p ;
                    this->post_prob_rowsum[i] += p ;
                }
            }
        }
        // normalize
        for(size_t n_class=0; n_class<this->n_class; n_class++)
        {   for(size_t n_shift=0; n_shift<this->n_shift; n_shift++)
            {   for(size_t n_flip=0; n_flip<this->n_flip; n_flip++)
                {
                    double p = std::max(this->post_prob(i,n_class,n_shift,n_flip) /
                                        this->post_prob_rowsum[i],
                                        ReadLayer::p_min) ;
                    // double p = this->post_prob(i,n_class,n_shift,n_flip) /
                    //            this->post_prob_rowsum[i] ;
                    this->post_prob(i,n_class,n_shift,n_flip) = p ;
                    colsums[n_class] += p ;
                }
            }
        }
    }
    post_prob_colsum.set_value(colsums) ;
}

void EMJoint::update_models()
{   // read data and models
    for(auto& ptr : this->read_layers)
    {   ptr->update_model(this->post_prob,
                          this->post_prob_colsum,
                          this->threads) ;
    }
    // sequence data and models
    if(this->seq_layer != nullptr)
    {   this->seq_layer->update_model(this->post_prob,
                                      this->threads) ;
    }
    // consensus sequence data and models
    if(this->consseq_layer != nullptr)
    {   this->consseq_layer->update_model(this->post_prob,
                                          this->threads) ;
    }
}
