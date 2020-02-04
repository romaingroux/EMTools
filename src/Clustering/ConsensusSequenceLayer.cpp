#include <ConsensusSequenceLayer.hpp>

#include <vector>
#include <future>      // std::promise, std::future
#include <stdexcept>   // std::invalid_argument
#include <cmath>       // log()
#include <functional>  // std::bind(), std::ref()

#include <ThreadPool.hpp>
#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <dna_utility.hpp>  // dna::base_composition()


double ConsensusSequenceLayer::score_consensus_subseq(const Matrix3D<double>& cons_seq,
                                                      size_t row,
                                                      size_t start,
                                                      const Matrix2D<double>& model)
{   size_t ncol = cons_seq.get_dim()[1] ;
    size_t dim3 = cons_seq.get_dim()[2] ;

    if(start > ncol - model.get_nrow())
    {   char msg[4096] ;
        sprintf(msg, "Error! given start (%zu) is too high. Max value is %zu",
                start, ncol - model.get_nrow()) ;
        throw std::invalid_argument(msg) ;
    }
    else if(model.get_nrow() > ncol)
    {   char msg[4096] ;
        sprintf(msg, "Error! given model is longer than sequences (%zu / %zu)",
                model.get_nrow(), ncol) ;
        throw std::invalid_argument(msg) ;
    }
    else if(model.get_ncol() != 4)
    {   char msg[4096] ;
        sprintf(msg, "Error! given model 2nd dimension is not 4 (%zu)",
                model.get_ncol()) ;
        throw std::invalid_argument(msg) ;
    }
    else if(dim3 != 4)
    {   char msg[4096] ;
        sprintf(msg, "Error! given data 3rd dimension is not 4 (%zu)",
                dim3) ;
        throw std::invalid_argument(msg) ;
    }

    size_t from = start ;
    size_t to   = from + model.get_nrow() ; // will score [from, to)


    double ll = 0 ;

    for(size_t j_seq=start, j_mod=0; j_seq<to; j_seq++, j_mod++)
    {   double sum = 0. ;
        for(size_t k=0; k<dim3; k++)
        {   sum += (cons_seq(row, j_seq, k) * model(j_mod, k)) ; }
        ll += log(sum) ;
    }

    return ll ;
}


ConsensusSequenceLayer::ConsensusSequenceLayer(const Matrix3D<double>& data,
                                               size_t n_class,
                                               size_t n_shift,
                                               bool flip,
                                               bool bckg_class)
    : Data3DLayer(data, n_class, n_shift, flip,bckg_class)
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

ConsensusSequenceLayer::ConsensusSequenceLayer(Matrix3D<double>&& data,
                                               size_t n_class,
                                               size_t n_shift,
                                               bool flip,
                                               bool bckg_class)
    : Data3DLayer(std::move(data), n_class, n_shift, flip, bckg_class)
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

ConsensusSequenceLayer::ConsensusSequenceLayer(const Matrix3D<double>& data,
                                               const Matrix3D<double>& model,
                                               bool flip,
                                               bool bckg_class)
    : Data3DLayer(data, model,flip, bckg_class)
{}

ConsensusSequenceLayer::ConsensusSequenceLayer(Matrix3D<double>&& data,
                                               Matrix3D<double>&& model,
                                               bool flip,
                                               bool bckg_class)
    : Data3DLayer(std::move(data), std::move(model),flip, bckg_class)
{}

ConsensusSequenceLayer::~ConsensusSequenceLayer()
{ ; }

void ConsensusSequenceLayer::compute_loglikelihoods(Matrix4D<double>& loglikelihood,
                                                    vector_d& loglikelihood_max,
                                                    ThreadPool* threads) const
{
    // dimension checks
    this->check_loglikelihood_dim(loglikelihood) ;
    this->check_loglikelihood_max_dim(loglikelihood_max) ;

    // compute the prob model and the prob reverse-complement model
    std::vector<Matrix2D<double>> model(this->n_class,
                                        Matrix2D<double>(this->l_model,
                                                         this->n_category,
                                                         0.)) ;
    std::vector<Matrix2D<double>> model_rev = model ;
    for(size_t i=0; i<this->n_class; i++)
    {   for(size_t j=0; j<this->l_model; j++)
        {   for(size_t k=0; k<this->n_category; k++)
            {   // forward
                model[i](j,k) = this->model(i,j,k) ;
                // reverse
                model_rev[i](this->l_model-j-1,this->n_category-k-1)
                        = this->model(i,j,k) ;
            }
        }
    }

    // don't parallelize
    if(threads == nullptr)
    {   std::promise<bool> promise ;
        std::future<bool> future = promise.get_future() ;
        this->compute_loglikelihoods_routine(0, this->n_dim1,
                                             loglikelihood,
                                             loglikelihood_max,
                                             model,
                                             model_rev,
                                             promise) ;
        future.get() ;
    }
    // parallelize
    else
    {   size_t n_threads = threads->getNThread() ;
        // compute the slices on which each thread will work
        std::vector<std::pair<size_t,size_t>> slices =
                ThreadPool::split_range(0, this->n_dim1, n_threads) ;

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
                                std::bind(&ConsensusSequenceLayer::compute_loglikelihoods_routine,
                                          this,
                                          slice.first,
                                          slice.second,
                                          std::ref(loglikelihood),
                                          std::ref(loglikelihood_max),
                                          std::ref(model),
                                          std::ref(model_rev),
                                          std::ref(promises[i])))) ;
        }
        // wait until all threads are done working
        for(auto& future : futures)
        {   future.get() ; }
        // -------------------------- threads stop ---------------------------
    }
}

void ConsensusSequenceLayer::compute_loglikelihoods_routine(size_t from,
                                                            size_t to,
                                                            Matrix4D<double>& loglikelihood,
                                                            vector_d& loglikelihood_max,
                                                            const std::vector<Matrix2D<double>>& model,
                                                            const std::vector<Matrix2D<double>>& model_rev,
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
                {   double ll_fw = score_consensus_subseq(this->data, i, s, model[j]) ;
                    loglikelihood(i,j,s,flip_states::FORWARD) = ll_fw  ;
                    // keep track of max per row
                    if(ll_fw > loglikelihood_max[i])
                    {   loglikelihood_max[i] = ll_fw ; }

                }
                // reverse
                if(this->flip)
                {   double ll_rev = score_consensus_subseq(this->data, i, s, model_rev[j]) ;
                    loglikelihood(i,j,s,flip_states::REVERSE) = ll_rev  ;
                    // keep track of max per row
                    if(ll_rev > loglikelihood_max[i])
                    {   loglikelihood_max[i] = ll_rev ; }
                }
            }
        }
    }
    done.set_value(true) ;
}

void ConsensusSequenceLayer::update_model(const Matrix4D<double>& posterior_prob,
                                          ThreadPool* threads)
{   // don't parallelize
    if(threads == nullptr)
    {   std::promise<Matrix3D<double>> promise ;
        std::future<Matrix3D<double>>  future = promise.get_future() ;
        this->update_model_routine(0,
                                   this->n_dim1,
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
                ThreadPool::split_range(0, this->n_dim1, n_threads) ;

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
                                std::bind(&ConsensusSequenceLayer::update_model_routine,
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
                        std::max(this->model(i,j,k), ConsensusSequenceLayer::p_min) ;
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

void ConsensusSequenceLayer::update_model_routine(size_t from,
                                                  size_t to,
                                                  const Matrix4D<double>& posterior_prob,
                                                  std::promise<Matrix3D<double>>& promise) const
{   // dimension checks
    this->check_posterior_prob_dim(posterior_prob) ;

    Matrix3D<double> model(this->n_class,
                           this->l_model,
                           this->n_category,
                           0.) ;

    size_t n_class_to_update = this->n_class - this->bckg_class ;

    for(size_t k=0; k < n_class_to_update; k++)
    {   for(size_t s=0; s<this->n_shift; s++)
        {   for(size_t j=0; j<this->l_model; j++)
            {   // base prob on fw and rv strand
                vector_d base_prob_fw(this->n_category, 0.) ;
                vector_d base_prob_rv(this->n_category, 0.) ;
                for(size_t i=from; i<to; i++)
                {   for(size_t base=0, base_rv=3; base<this->n_dim3; base++, base_rv--)
                    {   // --------------- forward ---------------
                        {   base_prob_fw[base] +=
                                    this->data(i,s+j,base) *
                                    posterior_prob(i,k,s,ConsensusSequenceLayer::FORWARD) ;
                        }
                        // --------------- reverse ---------------
                        if(this->flip)
                        {
                            base_prob_rv[base_rv] +=
                                    this->data(i,s+j,base) *
                                    posterior_prob(i,k,s,ConsensusSequenceLayer::REVERSE) ;
                        }
                    }
                }
                // update this position of the model
                for(size_t i=0,i_rev=base_prob_fw.size()-1; i<base_prob_fw.size(); i++,i_rev--)
                {   // --------------- forward ---------------
                    {   model(k,j,i) += base_prob_fw[i] ; }
                    // --------------- reverse ---------------
                    if(this->flip)
                    {   model(k,this->l_model-j-1,i) += base_prob_rv[i] ; }
                }
            }
        }
    }
    promise.set_value(model) ;
}

void ConsensusSequenceLayer::create_bckg_class()
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
