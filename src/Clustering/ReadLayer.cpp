#include <ReadLayer.hpp>

#include <stdexcept>         // std::invalid_argument
#include <limits>            // numeric_limits
#include <cmath>             // log(), exp(), pow()
#include <vector>
#include <future>            // std::promise, std::future
#include <utility>           // std::pair, std::move()
#include <functional>        // std::bind(), std::ref()

#include <Statistics.hpp>    // poisson_pmf()
#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

#include <iostream>

typedef std::vector<double> vector_d ;

ReadLayer::ReadLayer(const Matrix2D<int>& data,
                     size_t n_class,
                     size_t n_shift,
                     bool flip,
                     bool bckg_class,
                     ThreadPool* threads)
    : Data2DLayer(data, n_class, n_shift, flip, bckg_class),
      window_means(n_row, n_shift, 0.)
{   this->n_category = 1 ;
    // initialise the empty model
    this->model = Matrix3D<double>(this->n_class,
                                   this->l_model,
                                   this->n_category,
                                   0) ;
    // background class
    if(this->bckg_class)
    {   this->create_bckg_class() ; }

    // compute window means
    this->compute_window_means(threads) ;
}

ReadLayer::ReadLayer(Matrix2D<int>&& data,
                     size_t n_class,
                     size_t n_shift,
                     bool flip,
                     bool bckg_class,
                     ThreadPool* threads)
    : Data2DLayer(std::move(data), n_class, n_shift, flip, bckg_class),
      window_means(n_row, n_shift, 0.)
{   this->n_category = 1 ;
    // initialise the empty model
    this->model = Matrix3D<double>(this->n_class,
                                   this->l_model,
                                   this->n_category,
                                   0) ;

    // background class
    if(this->bckg_class)
    {   this->create_bckg_class() ; }

    // compute window means
    this->compute_window_means(threads) ;
}

ReadLayer::ReadLayer(const Matrix2D<int>& data,
                     const Matrix3D<double>& model,
                     bool flip,
                     bool bckg_class,
                     ThreadPool* threads)
    : Data2DLayer(data, model, flip, bckg_class),
      window_means(n_row, n_shift, 0.)
{   // check that the model only has one category
    if(this->n_category > 1)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! model is expected to have length 1 on "
                "3rd dimension, not %zu",
                this->n_category) ;
        throw std::invalid_argument(msg) ;
    }
    // compute window means
    this->compute_window_means(threads) ;
}

ReadLayer::ReadLayer(Matrix2D<int>&& data,
                     Matrix3D<double>&& model,
                     bool flip,
                     bool bckg_class,
                     ThreadPool* threads)
    : Data2DLayer(std::move(data), std::move(model), flip, bckg_class),
      window_means(n_row, n_shift, 0.)
{   // check that the model only has one category
    if(this->n_category > 1)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! model is expected to have length 1 on "
                "3rd dimension, not %zu",
                this->n_category) ;
        throw std::invalid_argument(msg) ;
    }
    // compute window means
    this->compute_window_means(threads) ;
}


ReadLayer::~ReadLayer()
{}

void ReadLayer::compute_loglikelihoods(Matrix4D<double>& loglikelihood,
                                       vector_d& loglikelihood_max,
                                       ThreadPool* threads) const
{   // dimension checks
    this->check_loglikelihood_dim(loglikelihood) ;
    this->check_loglikelihood_max_dim(loglikelihood_max) ;

    // don't parallelize
    if(threads == nullptr)
    {   std::promise<bool> promise ;
        std::future<bool> future = promise.get_future() ;
        this->compute_loglikelihoods_routine(0,
                                             this->n_row,
                                             std::ref(loglikelihood),
                                             std::ref(loglikelihood_max),
                                             promise) ;
        future.get() ;
    }
    // parallelize
    else
    {   size_t n_threads = threads->getNThread() ;
        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_row, n_threads) ;

        // get promises and futures
        // the function run by the threads will simply fill the promise with
        // "true" to indicate that they are done
        std::vector<std::promise<bool>> promises(n_threads) ;
        std::vector<std::future<bool>>  futures(n_threads) ;
        for(size_t i=0; i<n_threads; i++)
        {   futures[i] = promises[i].get_future() ; }
        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t i=0; i<n_threads; i++)
        {   auto slice = slices[i] ;
            threads->addJob(std::move(
                                std::bind(&ReadLayer::compute_loglikelihoods_routine,
                                          this,
                                          slice.first,
                                          slice.second,
                                          std::ref(loglikelihood),
                                          std::ref(loglikelihood_max),
                                          std::ref(promises[i])))) ;
        }
        // wait until all threads are done working
        for(auto& future : futures)
        {   future.get() ; }
        // -------------------------- threads stop ---------------------------
    }
}


void ReadLayer::compute_loglikelihoods_routine(size_t from,
                                               size_t to,
                                               Matrix4D<double>& loglikelihood,
                                               vector_d& loglikelihood_max,
                                               std::promise<bool>& done) const
{
    // normalize the models
    Matrix3D<double> model_norm = this->model ;
    for(size_t i=0; i<this->n_class; i++)
    {   double mean = 0. ;
        for(size_t j=0; j<this->l_model; j++)
        {   mean += model_norm(i,j,0) ; }
        mean /= this->l_model ;
        for(size_t j=0; j<this->l_model; j++)
        {   model_norm(i,j,0) /= mean ; }
    }

    // compute log likelihood
    for(size_t i=from; i<to; i++)
    {
        // set max to min possible value
        loglikelihood_max[i] = std::numeric_limits<double>::lowest() ;

        for(size_t j=0; j<this->n_class; j++)
        {   for(size_t s_fw=0, s_rev=this->n_shift-1;
                s_fw<this->n_shift; s_fw++, s_rev--)
            {   // slice is [from_fw,to)
                //    from_dat_fw             to_dat_fw    [from_dat_fw, to_dat_fw]
                // fw  |---------->>>----------|
                //     ----------------------------------> data
                // rev           |----------<<<----------| [from_dat_rev, to_dat_rev]
                //                                         to_dat_rev can be -1 -> int
                //            to_dat_rev             from_dat_rev

                // log likelihood
                double ll_fw  = 0. ;
                double ll_rev = 0. ;
                // --------------- forward ---------------
                size_t from_dat_fw   = s_fw ;
                size_t to_dat_fw     = from_dat_fw + this->l_model - 1 ;
                // --------------- reverse ---------------
                size_t from_dat_rev = this->n_col - 1 - s_fw ;
                // size_t to_dat_rev   = from_dat_rev - (this->l_model - 1) ;

                for(size_t j_dat_fw=from_dat_fw,j_ref_fw=0, j_dat_rev=from_dat_rev;
                    j_dat_fw<to_dat_fw;
                    j_dat_fw++, j_ref_fw++, j_dat_rev--)
                {
                    double ll ;
                    // --------------- forward ---------------
                    ll = log(poisson_pmf(this->data(i,j_dat_fw),
                                         model_norm(j,j_ref_fw,0)*
                                         this->window_means(i,s_fw))) ;
                    // p(A|B) may be really unlikely -> rounded to 0 -> log(0) = -inf
                    // prevent this by applying this
                    if(isinf(ll))
                    {   ll = ReadLayer::p_min_log ; }
                    ll_fw += ll ;
                    // --------------- reverse ---------------
                    if(this->flip)
                    {   ll = log(poisson_pmf(this->data(i,j_dat_rev),
                                             model_norm(j,j_ref_fw,0)*
                                             this->window_means(i,s_rev))) ;
                        // p(A|B) may be really unlikely -> rounded to 0 -> log(0) = -inf
                        // prevent this by applying this
                        if(isinf(ll))
                        {   ll = ReadLayer::p_min_log ; }
                        ll_rev += ll ;
                    }
                }
                loglikelihood(i,j,from_dat_fw,flip_states::FORWARD) = ll_fw  ;
                // keep track of the max per row
                if(ll_fw > loglikelihood_max[i])
                {   loglikelihood_max[i] = ll_fw ; }

                if(this->flip)
                {   loglikelihood(i,j,from_dat_fw,flip_states::REVERSE) = ll_rev ;
                    // keep track of the max per row
                    if(ll_rev > loglikelihood_max[i])
                    {   loglikelihood_max[i] = ll_rev ; }
                }
            }
        }
    }
    done.set_value(true) ;
}


void ReadLayer::update_model(const Matrix4D<double>& posterior_prob,
                             ThreadPool* threads)
{
    // computing sum over the columns (classes)
    size_t n_row   = posterior_prob.get_dim()[0] ;
    size_t n_class = posterior_prob.get_dim()[1] ;
    size_t n_shift = posterior_prob.get_dim()[2] ;
    size_t n_flip  = posterior_prob.get_dim()[3] ;
    vector_d colsum(n_class, 0.) ;
    for(size_t i=0; i<n_row; i++)
    {   for(size_t j=0; j<n_class; j++)
        {   for(size_t k=0; k<n_shift; k++)
            {   for(size_t l=0; l<n_flip; l++)
                {   colsum[j] += posterior_prob(i,j,k,l) ; }
            }
        }
    }
    this->update_model(posterior_prob,
                       colsum,
                       threads) ;
}

void ReadLayer::update_model(const Matrix4D<double>& posterior_prob,
                             const vector_d& posterior_prob_colsum,
                             ThreadPool* threads)
{
    // don't parallelize
    if(threads == nullptr)
    {   std::promise<Matrix3D<double>> promise ;
        std::future<Matrix3D<double>>  future = promise.get_future() ;
        this->update_model_routine(0,
                                   this->n_row,
                                   posterior_prob,
                                   posterior_prob_colsum,
                                   promise) ;
        // this->model = future.get() ;
        auto model = future.get() ;
        size_t n_class_to_update = this->n_class - this->bckg_class ;
        for(size_t i=0; i<n_class_to_update; i++)
        {   for(size_t j=0; j<this->l_model; j++)
            {   for(size_t k=0; k<this->n_category; k++)
                {   this->model(i,j,k) = model(i,j,k) ; }
            }
        }
    }
    // parallelize
    else
    {   size_t n_threads = threads->getNThread() ;
        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_row, n_threads) ;

        // get promises and futures
        // the function run by the threads will simply fill the promise with
        // "true" to indicate that they are done
        std::vector<std::promise<Matrix3D<double>>> promises(n_threads) ;
        std::vector<std::future<Matrix3D<double>>>  futures(n_threads) ;
        for(size_t i=0; i<n_threads; i++)
        {   futures[i] = promises[i].get_future() ; }
        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t i=0; i<n_threads; i++)
        {   auto slice = slices[i] ;
            threads->addJob(std::move(
                                std::bind(&ReadLayer::update_model_routine,
                                          this,
                                          slice.first,
                                          slice.second,
                                          std::ref(posterior_prob),
                                          std::ref(posterior_prob_colsum),
                                          std::ref(promises[i])))) ;
        }
        // reinitialise the model
        /*
        this->model = Matrix3D<double>(this->n_class,
                                       this->l_model,
                                       this->n_category,
                                       0.) ;
        */
        size_t n_class_to_update = this->n_class - this->bckg_class ;
        for(size_t i=0; i<n_class_to_update; i++)
        {   for(size_t j=0; j<this->l_model; j++)
            {   for(size_t k=0; k<this->n_category; k++)
                {   this->model(i,j,k) = 0. ; }
            }
        }
        // wait until all threads are done working
        // and update the mode
        for(auto& future : futures)
        {   Matrix3D<double> model_part = future.get() ;
            // for(size_t i=0; i<this->n_class; i++)
            for(size_t i=0; i<n_class_to_update; i++)
            {   for(size_t j=0; j<this->l_model; j++)
                {   for(size_t k=0; k<this->n_category; k++)
                    {   this->model(i,j,k) +=
                            model_part(i,j,k) ;
                    }
                }
            }
        }
        // -------------------------- threads stop ---------------------------
    }
    // avoid 0's in the model to ensure that pmf_poisson() never
    // return 0
    for(size_t i=0; i<this->n_class; i++)
    {   for(size_t j=0; j<this->l_model; j++)
        {   for(size_t k=0; k<this->n_category; k++)
            {   this->model(i,j,k) =
                    std::max(this->model(i,j,k), ReadLayer::p_min) ;
            }
        }
    }
}

void ReadLayer::update_model_routine(size_t from,
                                     size_t to,
                                     const Matrix4D<double>& posterior_prob,
                                     const vector_d& posterior_prob_colsum,
                                     std::promise<Matrix3D<double>>& promise) const
{
    // dimension checks
    this->check_posterior_prob_dim(posterior_prob) ;
    this->check_posterior_prob_colsum_dim(posterior_prob_colsum) ;

    // partial model
    Matrix3D<double> model(this->n_class,
                           this->l_model,
                           this->n_category,
                           0.) ;

    size_t n_class_to_update = this->n_class - this->bckg_class ;

    for(size_t n_class=0; n_class < n_class_to_update; n_class++)
    {   for(size_t i=from; i<to; i++)
        {   for(size_t n_shift=0; n_shift<this->n_shift; n_shift++)
            {   // --------------- forward ---------------
                int from_dat_fw = n_shift ;
                int to_dat_fw   = from_dat_fw + this->l_model - 1 ;
                for(int j_dat_fw=from_dat_fw, j_ref_fw=0;
                    j_dat_fw<=to_dat_fw; j_dat_fw++, j_ref_fw++)
                {   model(n_class,j_ref_fw,0) +=
                            (posterior_prob(i,n_class,n_shift,flip_states::FORWARD) *
                            this->data(i,j_dat_fw)) /
                            posterior_prob_colsum[n_class] ;
                }
                // --------------- reverse ---------------
                if(this->flip)
                {   int from_dat_rev = this->n_col - 1 - n_shift ;
                    int to_dat_rev   = from_dat_rev - (this->l_model - 1) ;
                    for(int j_dat_rev=from_dat_rev, j_ref_fw=0;
                        j_dat_rev >= to_dat_rev; j_dat_rev--, j_ref_fw++)
                    {   model(n_class,j_ref_fw,0) +=
                                (posterior_prob(i,n_class,n_shift,flip_states::REVERSE) *
                                this->data(i,j_dat_rev)) /
                                posterior_prob_colsum[n_class] ;
                    }
                }
            }
        }
    }
    promise.set_value(model) ;
}

void ReadLayer::compute_window_means(ThreadPool* threads)
{   // don't parallelize
    if(threads == nullptr)
    {   std::promise<bool> promise ;
        std::future<bool> future = promise.get_future() ;
        this->compute_window_means_routine(0, this->n_row, promise) ;
        future.get() ;
    }
    // parallelize
    else
    {   size_t n_threads = threads->getNThread() ;
        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_row, n_threads) ;

        // get promises and futures
        // the function run by the threads will simply fill the promise with
        // "true" to indicate that they are done
        std::vector<std::promise<bool>> promises(n_threads) ;
        std::vector<std::future<bool>>  futures(n_threads) ;
        for(size_t i=0; i<n_threads; i++)
        {   futures[i] = promises[i].get_future() ; }
        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t i=0; i<n_threads; i++)
        {   auto slice = slices[i] ;
            threads->addJob(std::move(
                                std::bind(&ReadLayer::compute_window_means_routine,
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

void ReadLayer::compute_window_means_routine(size_t from,
                                             size_t to,
                                             std::promise<bool>& done)
{   double l_window = double(this->l_model) ;
    for(size_t i=from; i<to; i++)
    {   for(size_t from=0; from<this->n_shift; from++)
        {   double sum = 0. ;
            // slice is [from,to)
            size_t to = from + this->l_model ;
            for(size_t j=from; j<to; j++)
            {   sum += this->data(i,j) ;}
            this->window_means(i,from) = sum / l_window ;
        }
    }
    done.set_value(true) ;
}

void ReadLayer::check_posterior_prob_colsum_dim(const vector_d& posterior_prob_colsum) const
{   if(posterior_prob_colsum.size() != this->n_class)
    {   char msg[4096] ;
        sprintf(msg,
                "Error! posterior_class_prob matrix size is not "
                "equal to model class number : %zu / %zu",
                posterior_prob_colsum.size(), this->n_class) ;
        throw std::invalid_argument(msg) ;
    }
}

void ReadLayer::create_bckg_class()
{
    // mean count
    double mean = 0. ;
    double n = (double)this->data.get_nrow() *
               (double)this->data.get_ncol() ;
    for(size_t i=0; i<this->data.get_nrow(); i++)
    {   for(size_t j=0; j<this->data.get_ncol(); j++)
        {   mean += ((double)this->data(i,j)) / n ; }
    }
    // set last class as background class
    for(size_t j=0; j<this->l_model; j++)
    {   this->model(this->n_class-1,j,0) = mean ; }
}

