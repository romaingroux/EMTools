#include <SequenceLayer.hpp>

#include <stdexcept>         // std::invalid_argument
#include <limits>            // numeric_limits
#include <cmath>             // log(), pow()
#include <vector>
#include <algorithm>         // std::max_element()
#include <future>            // std::promise, std::future

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <dna_utility.hpp>  // dna::base_composition()


double SequenceLayer::score_subseq(const Matrix2D<int>& seq,
                                   size_t row,
                                   size_t start,
                                   const Matrix2D<double>& model_log)
{
    if(start > seq.get_ncol() - model_log.get_nrow())
    {   char msg[4096] ;
        sprintf(msg, "Error! given start (%zu) is too high. Max value is %zu",
                start, seq.get_ncol() - model_log.get_nrow()) ;
        throw std::invalid_argument(msg) ;
    }
    else if(model_log.get_nrow() > seq.get_ncol())
    {   char msg[4096] ;
        sprintf(msg, "Error! given model is longer than sequences (%zu / %zu)",
                model_log.get_nrow(), seq.get_ncol()) ;
        throw std::invalid_argument(msg) ;
    }
    else if(model_log.get_ncol() != 4)
    {   char msg[4096] ;
        sprintf(msg, "Error! given model 2nd dimension is not 4 (%zu)",
                model_log.get_ncol()) ;
        throw std::invalid_argument(msg) ;
    }

    size_t from = start ;
    size_t to   = from + model_log.get_nrow() ; // will score [from, to)

    int n_code = dna::char_to_int('N') ;
    double ll = 0 ;
    for(size_t i=from, j=0; i<to; i++, j++)
    {   int base = seq(row,i) ;
        // N char -> get max score
        if(base == n_code)
        {   std::vector<double> row = model_log.get_row(j) ;
            ll += *(std::max_element(std::begin(row),
                                     std::end(row)))  ;
        }
        // A,C,G,T -> get its score
        else
        {   ll += model_log(j,base) ; }
    }
    return ll ;
}

SequenceLayer::SequenceLayer(const Matrix2D<int>& data,
                             size_t n_class,
                             size_t n_shift,
                             bool flip,
                             bool bckg_class)
    : Data2DLayer(data, n_class, n_shift, flip,bckg_class)
{
    this->n_category = 4 ;

    // initialise the empty model
    this->model = Matrix3D<double>(this->n_class,
                                   this->l_model,
                                   this->n_category,
                                   0.) ;
    // background class
    if(this->bckg_class)
    {   this->create_bckg_class() ; }
}

SequenceLayer::SequenceLayer(Matrix2D<int>&& data,
                             size_t n_class,
                             size_t n_shift,
                             bool flip,
                             bool bckg_class)
    : Data2DLayer(std::move(data), n_class, n_shift, flip, bckg_class)
{
    this->n_category = 4 ;

    // initialise the empty model
    this->model = Matrix3D<double>(this->n_class,
                                   this->l_model,
                                   this->n_category,
                                   0.) ;

    // background class
    if(this->bckg_class)
    {   this->create_bckg_class() ; }
}

SequenceLayer::SequenceLayer(const Matrix2D<int>& data,
                             const Matrix3D<double>& model,
                             bool flip,
                             bool bckg_class)
    : Data2DLayer(data, model,flip, bckg_class)
{}

SequenceLayer::SequenceLayer(Matrix2D<int>&& data,
                             Matrix3D<double>&& model,
                             bool flip,
                             bool bckg_class)
    : Data2DLayer(std::move(data), std::move(model),flip, bckg_class)
{}

SequenceLayer::~SequenceLayer()
{ ; }

void SequenceLayer::compute_loglikelihoods(Matrix4D<double>& loglikelihood,
                                           vector_d& loglikelihood_max,
                                           ThreadPool* threads) const
{
    // dimension checks
    this->check_loglikelihood_dim(loglikelihood) ;
    this->check_loglikelihood_max_dim(loglikelihood_max) ;

    // compute the log prob model and the log prob reverse-complement model
    std::vector<Matrix2D<double>> model_log(this->n_class,
                                            Matrix2D<double>(this->l_model,
                                                             this->n_category,
                                                             0.)) ;
    std::vector<Matrix2D<double>> model_log_rev = model_log ;
    for(size_t i=0; i<this->n_class; i++)
    {   for(size_t j=0; j<this->l_model; j++)
        {   for(size_t k=0; k<this->n_category; k++)
            {   // forward
                model_log[i](j,k) = log(this->model(i,j,k)) ;
                // reverse
                model_log_rev[i](this->l_model-j-1,this->n_category-k-1)
                        = log(this->model(i,j,k)) ;
            }
        }
    }

    // don't parallelize
    if(threads == nullptr)
    {   std::promise<bool> promise ;
        std::future<bool> future = promise.get_future() ;
        this->compute_loglikelihoods_routine(0, this->n_row,
                                             loglikelihood,
                                             loglikelihood_max,
                                             model_log,
                                             model_log_rev,
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
                                std::bind(&SequenceLayer::compute_loglikelihoods_routine,
                                          this,
                                          slice.first,
                                          slice.second,
                                          std::ref(loglikelihood),
                                          std::ref(loglikelihood_max),
                                          std::ref(model_log),
                                          std::ref(model_log_rev),
                                          std::ref(promises[i])))) ;
        }
        // wait until all threads are done working
        for(auto& future : futures)
        {   future.get() ; }
        // -------------------------- threads stop ---------------------------
    }
}

void SequenceLayer::compute_loglikelihoods_routine(size_t from,
                                                   size_t to,
                                                   Matrix4D<double>& loglikelihood,
                                                   vector_d& loglikelihood_max,
                                                   const std::vector<Matrix2D<double>>& model_log,
                                                   const std::vector<Matrix2D<double>>& model_log_rev,
                                                   std::promise<bool>& done) const
{
    // compute log likelihood
    for(size_t i=from; i<to; i++)
    {
        // set max to min possible value
        loglikelihood_max[i] = std::numeric_limits<double>::lowest() ;

        for(size_t j=0; j<this->n_class; j++)
        {
            for(size_t s=0; s<this->n_shift; s++)
            {   // forward strand
                {   double ll_fw = score_subseq(this->data, i, s, model_log[j]) ;
                    // loglikelihood(i,j,s,flip_states::FORWARD) = ll_fw ;
                    // apply lower bound here -> log(1e-100)
                    loglikelihood(i,j,s,flip_states::FORWARD) = ll_fw  ;
                    // keep track of max per row
                    if(ll_fw > loglikelihood_max[i])
                    {   loglikelihood_max[i] = ll_fw ; }

                }
                // reverse
                if(this->flip)
                {   double ll_rev = score_subseq(this->data, i, s, model_log_rev[j]) ;
                    // loglikelihood(i,j,s,flip_states::REVERSE) = ll_rev ;
                    // apply lower bound here -> log(1e-100)
                    loglikelihood(i,j,s,flip_states::REVERSE) = ll_rev ;
                    // keep track of max per row
                    if(ll_rev > loglikelihood_max[i])
                    {   loglikelihood_max[i] = ll_rev ; }
                }
            }
        }
    }
    done.set_value(true) ;
}


void SequenceLayer::update_model(const Matrix4D<double>& posterior_prob,
                                 ThreadPool* threads)
{   // don't parallelize
    if(threads == nullptr)
    {   std::promise<Matrix3D<double>> promise ;
        std::future<Matrix3D<double>>  future = promise.get_future() ;
        this->update_model_routine(0,
                                   this->n_row,
                                   posterior_prob,
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
                                std::bind(&SequenceLayer::update_model_routine,
                                          this,
                                          slice.first,
                                          slice.second,
                                          std::ref(posterior_prob),
                                          std::ref(promises[i])))) ;
        }
        // reinitialise the model
        size_t n_class_to_update = this->n_class - this->bckg_class ;
        for(size_t i=0; i<n_class_to_update; i++)
        {   for(size_t j=0; j<this->l_model; j++)
            {   for(size_t k=0; k<this->n_category; k++)
                {   this->model(i,j,k) = 0. ; }
            }
        }
        // wait until all threads are done working
        // and update the model
        for(auto& future : futures)
        {   Matrix3D<double> model_part = future.get() ;
            // for(size_t i=0; i<this->n_class; i++)
            for(size_t i=0; i<n_class_to_update; i++)
            {   for(size_t j=0; j<this->l_model; j++)
                {   for(size_t k=0; k<this->n_category; k++)
                    {   this->model(i,j,k) += model_part(i,j,k) ; }
                }
            }
        }
        // -------------------------- threads stop ---------------------------
    }
    // make sure to have no 0 values
    for(size_t i=0; i<this->n_class; i++)
    {   for(size_t j=0; j<this->l_model; j++)
        {   for(size_t k=0; k<this->n_category; k++)
            {   this->model(i,j,k) =
                        std::max(this->model(i,j,k), SequenceLayer::p_min) ;
            }
        }
    }
    // normalize to get probs
    for(size_t i=0; i<this->n_class; i++)
    {   for(size_t j=0; j<this->l_model; j++)
        {   double sum = 0. ;
            for(size_t k=0; k<this->n_category; k++)
            {   sum += this->model(i,j,k) ; }
            for(size_t k=0; k<this->n_category; k++)
            {   double p = this->model(i,j,k) / sum ;
                this->model(i,j,k) = p ;
            }
        }
    }
}


void SequenceLayer::update_model_routine(size_t from,
                                         size_t to,
                                         const Matrix4D<double>& posterior_prob,
                                         std::promise<Matrix3D<double>>& promise) const
{   // dimension checks
    this->check_posterior_prob_dim(posterior_prob) ;

    Matrix3D<double> model(this->n_class,
                           this->l_model,
                           this->n_category,
                           0.) ;

    // the int code of A, C, G, T, N
    static int a_code = dna::char_to_int('A') ;
    static int c_code = dna::char_to_int('C') ;
    static int g_code = dna::char_to_int('G') ;
    static int t_code = dna::char_to_int('T') ;
    static int n_code = dna::char_to_int('N') ;
    // the int code of the reverse complement of A, C, G, T
    static int a_code_r = dna::char_to_int('A', true) ;
    static int c_code_r = dna::char_to_int('C', true) ;
    static int g_code_r = dna::char_to_int('G', true) ;
    static int t_code_r = dna::char_to_int('T', true) ;

    size_t n_class_to_update = this->n_class - this->bckg_class ;

    for(size_t k=0; k<n_class_to_update; k++)
    {   for(size_t s=0; s<this->n_shift; s++)
        {   // for(size_t j=0; j<this->l_model; j++)
            for(size_t j=0, j_rv=this->l_model-1; j<this->l_model; j++, j_rv--)
            {   // base prob on fw and rv strand
                vector_d base_prob_fw(this->n_category, 0.) ;
                vector_d base_prob_rv(this->n_category, 0.) ;
                for(size_t i=from; i<to; i++)
                {   int base     = this->data(i,s+j) ;
                    int base_rv = this->n_category - base - 1 ; // complement
                    // int base_rv = this->n_category - this->data(i,s+j_rv) - 1 ;
                    // std::cerr << k << " "
                    //           << s << " "
                    //           << j    << "/" << j_rv << " "
                    //           << base << "/" << base_rv
                    //           << std::endl ;
                    // N
                    if(base == n_code)
                    {   // --------------- forward ---------------
                        {   base_prob_fw[a_code] +=
                                    posterior_prob(i,k,s,SequenceLayer::FORWARD) ;
                            base_prob_fw[c_code] +=
                                    posterior_prob(i,k,s,SequenceLayer::FORWARD) ;
                            base_prob_fw[g_code] +=
                                    posterior_prob(i,k,s,SequenceLayer::FORWARD) ;
                            base_prob_fw[t_code] +=
                                    posterior_prob(i,k,s,SequenceLayer::FORWARD) ;
                        }
                        // --------------- reverse ---------------
                        if(this->flip)
                        {   base_prob_rv[a_code_r] +=
                                    posterior_prob(i,k,s,SequenceLayer::REVERSE) ;
                            base_prob_rv[c_code_r] +=
                                    posterior_prob(i,k,s,SequenceLayer::REVERSE) ;
                            base_prob_rv[g_code_r] +=
                                    posterior_prob(i,k,s,SequenceLayer::REVERSE) ;
                            base_prob_rv[t_code_r] +=
                                    posterior_prob(i,k,s,SequenceLayer::REVERSE) ;
                        }
                    }
                    // A, C, G, T
                    else
                    {   // --------------- forward ---------------
                        {   base_prob_fw[base] +=
                                    posterior_prob(i,k,s,SequenceLayer::FORWARD) ;
                        }
                        // --------------- reverse ---------------
                        if(this->flip)
                        {   base_prob_rv[base_rv] +=
                                    posterior_prob(i,k,s,SequenceLayer::REVERSE) ;
                        }
                    }
                }
                // update this position of the model
                for(size_t i=0; i<base_prob_fw.size(); i++)
                {   // --------------- forward ---------------
                    {   model(k,j,i) += base_prob_fw[i] ; }
                    // --------------- reverse ---------------
                    if(this->flip)
                    {   model(k,this->l_model-j-1,i) += base_prob_rv[i] ;
                        /* model(k,j_rv,i) += base_prob_rv[i] ; */

                    }
                }
            }
        }
    }
    promise.set_value(model) ;
}

void SequenceLayer::create_bckg_class()
{   // sequence composition
    std::vector<double> base_comp =
                      dna::base_composition(this->data,
                                            this->flip) ;
    // set last class as background class
    for(size_t i=0; i<this->n_category; i++)
    {   for(size_t j=0; j<this->l_model; j++)
        {   this->model(this->n_class-1,j,i) = base_comp[i] ; }
    }
}
