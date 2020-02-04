#include <EMRead.hpp>

#include <string>
#include <vector>
#include <future>                    // std::promise, std::future
#include <utility>                   // std::pair, std::move()
#include <functional>                // std::bind(), std::ref()
#include <cmath>                     // exp()

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <ReadLayer.hpp>             // ReadLayer
#include <RandomNumberGenerator.hpp> // getRandomNumberGenerator()
#include <ConsoleProgressBar.hpp>    // ConsoleProgressBar
#include <ThreadPool.hpp>            // ThreadPool


EMRead::EMRead(const Matrix2D<int>& read_matrix,
               size_t n_class,
               size_t n_iter,
               size_t n_shift,
               bool flip,
               bool bckg_class,
               const std::string& seed,
               size_t n_threads)
    : EMBase(read_matrix.get_nrow(),
             read_matrix.get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      loglikelihood_max(n_row, 0.),
      read_layer(nullptr)
{   this->loglikelihood_max = vector_d(n_row, 0.) ;

    // initialise post prob randomly
    this->set_post_prob_random(seed) ;
    // data and models
    this->read_layer = new ReadLayer(read_matrix,
                                     this->n_class,
                                     this->n_shift,
                                     flip,
                                     bckg_class,
                                     this->threads) ;
    // intialise the models with the post prob
    this->read_layer->update_model(this->post_prob,
                                  this->threads) ;
}

EMRead::EMRead(Matrix2D<int>&& read_matrix,
               size_t n_class,
               size_t n_iter,
               size_t n_shift,
               bool flip,
               bool bckg_class,
               const std::string& seed,
               size_t n_threads)
    : EMBase(read_matrix.get_nrow(),
             read_matrix.get_ncol(),
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      loglikelihood_max(n_row, 0.),
      read_layer(nullptr)
{   this->loglikelihood_max = vector_d(n_row, 0.) ;

    // initialise post prob randomly
    this->set_post_prob_random(seed) ;

    // data and models
    this->read_layer = new ReadLayer(std::move(read_matrix),
                                     this->n_class,
                                     this->n_shift,
                                     flip,
                                     bckg_class,
                                     this->threads) ;
    // intialise the models with the post prob
    this->read_layer->update_model(this->post_prob,
                                  this->threads) ;
}

EMRead::~EMRead()
{   if(this->read_layer!= nullptr)
    {   delete this->read_layer ;
        this->read_layer = nullptr ;
    }
    if(this->threads != nullptr)
    {   this->threads->join() ;
        delete this->threads ;
        this->threads = nullptr ;
    }
}

Matrix3D<double> EMRead::get_read_models() const
{   return read_layer->get_model() ; }

EMRead::exit_codes EMRead::classify()
{
    size_t bar_update_n = this->n_iter ;
    ConsoleProgressBar bar(std::cerr, bar_update_n, 60, "classifying") ;

    // optimize the partition
    for(size_t n_iter=0; n_iter<this->n_iter; n_iter++)
    {
        // E-step
        this->compute_loglikelihood() ;
        // std::cerr << this->post_prob_rowsum << std::endl ;
        // std::cerr << this->post_prob_colsum << std::endl ;
        this->compute_post_prob() ;
        // M-step
        // std::cerr << this->post_prob_rowsum << std::endl ;
        // std::cerr << this->post_prob_colsum << std::endl ;
        this->compute_class_prob() ;
        this->update_models() ;
        this->center_post_state_prob() ;

        bar.update() ;
    }
    bar.update() ; std::cerr << std::endl ;
    return EMRead::exit_codes::ITER_MAX ;
}

void EMRead::compute_loglikelihood()
{   // compute the loglikelihood
    this->read_layer->compute_loglikelihoods(this->loglikelihood,
                                             this->loglikelihood_max,
                                             this->threads) ;

    // rescale the values
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
    {    size_t n_threads = this->threads->getNThread() ;

        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_row,n_threads) ;

        // get promises and futures
        std::vector<std::promise<bool>> promises(n_threads) ;
        std::vector<std::future<bool>>  futures(n_threads) ;
        for(size_t i=0; i<n_threads; i++)
        {   futures[i] = promises[i].get_future() ; }

        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t i=0; i<n_threads; i++)
        {   auto slice = slices[i] ;
            this->threads->addJob(std::move(
                                      std::bind(&EMRead::compute_loglikelihood_routine,
                                                this,
                                                slice.first,
                                                slice.second,
                                                std::ref(promises[i])))) ;
        }
        // wait until all threads are done working
        for(auto& future : futures)
        {   future.get() ; }
        // -------------------------- threads stop ---------------------------
    }
}

void EMRead::compute_loglikelihood_routine(size_t from,
                                           size_t to,
                                           std::promise<bool>& done)
{
    // rescale the values
    for(size_t i=from; i<to; i++)
    {   for(size_t j=0; j<this->n_class; j++)
        {   for(size_t k=0; k<this->n_shift; k++)
            {   for(size_t l=0; l<this->n_flip; l++)
                {   this->loglikelihood(i,j,k,l) =
                            std::max(this->loglikelihood(i,j,k,l) -
                                     this->loglikelihood_max[i],
                                     ReadLayer::p_min_log) ;
                }
            }
        }
    }
    done.set_value(true) ;
}

void EMRead::compute_post_prob()
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
                                      std::bind(&EMRead::compute_post_prob_routine,
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


void EMRead::compute_post_prob_routine(size_t from,
                                         size_t to,
                                         std::promise<vector_d>& post_prob_colsum)
{   vector_d colsums(this->n_class, 0.) ;

    // reset grand total
    // this->post_prob_tot = 0 ;
    // this->post_prob_colsum = vector_d(n_class, 0) ;

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
                    this->post_prob(i,n_class,n_shift,n_flip) = p ;
                    this->post_prob_rowsum[i] += p ;
                }
            }
        }

        // normalize
        for(size_t n_class=0; n_class<this->n_class; n_class++)
        {   for(size_t n_shift=0; n_shift<this->n_shift; n_shift++)
            {   for(size_t n_flip=0; n_flip<this->n_flip; n_flip++)
                {   // avoid p=0. by rounding errors
                    double p = std::max(this->post_prob(i,n_class,n_shift,n_flip) /
                                        this->post_prob_rowsum[i],
                                        ReadLayer::p_min) ;
                    this->post_prob(i,n_class,n_shift,n_flip) = p ;
                    colsums[n_class] += p ;
                }
            }
        }
    }
    post_prob_colsum.set_value(colsums) ;
}

void EMRead::update_models()
{   this->read_layer->update_model(this->post_prob,
                                   this->post_prob_colsum,
                                   this->threads) ;
}
