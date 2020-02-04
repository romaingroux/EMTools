#include <EMConsensusSequence.hpp>

#include <string>
#include <vector>
#include <future>                    // std::promise, std::future
#include <utility>                   // std::pair, std::move()
#include <functional>                // std::bind(), std::ref()

#include <ConsensusSequenceLayer.hpp>  // SequenceLayer
#include <RandomNumberGenerator.hpp>   // getRandomNumberGenerator()
#include <ConsoleProgressBar.hpp>      // ConsoleProgressBar
#include <ThreadPool.hpp>              // ThreadPool
#include <dna_utility.hpp>             // dna::base_composition()


EMConsensusSequence::EMConsensusSequence(const Matrix3D<double>& seq_matrix,
                                         size_t n_class,
                                         size_t n_iter,
                                         size_t n_shift,
                                         bool flip,
                                         bool bckg_class,
                                         const std::string& seed,
                                         size_t n_threads)
    : EMBase(seq_matrix.get_dim()[0],
             seq_matrix.get_dim()[1],
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      loglikelihood_max(n_row, 0.),
      cseq_layer(nullptr)
{
    this->loglikelihood_max = vector_d(n_row, 0.) ;

    // initialise post prob randomly
    // getRandomGenerator(seed) ;
    this->set_post_prob_random(seed) ;

    // data and models
    this->cseq_layer = new ConsensusSequenceLayer(seq_matrix,
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class) ;

    // intialise the models with the post prob
    this->cseq_layer->update_model(this->post_prob,
                                   this->threads) ;
}

EMConsensusSequence::EMConsensusSequence(Matrix3D<double>&& seq_matrix,
                                         size_t n_class,
                                         size_t n_iter,
                                         size_t n_shift,
                                         bool flip,
                                         bool bckg_class,
                                         const std::string& seed,
                                         size_t n_threads)
    : EMBase(seq_matrix.get_dim()[0],
             seq_matrix.get_dim()[1],
             n_class,
             n_iter,
             n_shift,
             flip,
             n_threads),
      loglikelihood_max(n_row, 0.),
      cseq_layer(nullptr)
{
    this->loglikelihood_max = vector_d(n_row, 0.) ;

    // initialise post prob randomly
    // getRandomGenerator(seed) ;
    this->set_post_prob_random(seed) ;

    // data and models
    this->cseq_layer = new ConsensusSequenceLayer(std::move(seq_matrix),
                                                  this->n_class,
                                                  this->n_shift,
                                                  this->flip,
                                                  bckg_class) ;

    // intialise the models with the post prob
    this->cseq_layer->update_model(this->post_prob,
                                   this->threads) ;
}

EMConsensusSequence::EMConsensusSequence(const Matrix3D<double>& seq_matrix,
                                         const Matrix3D<double>& motifs,
                                         size_t n_iter,
                                         bool flip,
                                         bool bckg_class,
                                         size_t n_threads)
    : EMBase(seq_matrix.get_dim()[0],
             seq_matrix.get_dim()[1],
             motifs.get_dim()[0],
             n_iter,
             seq_matrix.get_dim()[1] - motifs.get_dim()[1] + 1,
             flip,
             n_threads),
      loglikelihood_max(n_row, 0.),
      cseq_layer(nullptr)
{
    this->loglikelihood_max = vector_d(n_row, 0.) ;

    // data and models
    // background motif (if any) is the last of the given motifs
    this->cseq_layer = new ConsensusSequenceLayer(seq_matrix,
                                                  motifs,
                                                  this->flip,
                                                  bckg_class) ;

    // intialise the class prob uniformly
    this->set_state_prob_uniform() ;
}

EMConsensusSequence::EMConsensusSequence(Matrix3D<double>&& seq_matrix,
                                         Matrix3D<double>&& motifs,
                                         size_t n_iter,
                                         bool flip,
                                         bool bckg_class,
                                         size_t n_threads)
    : EMBase(seq_matrix.get_dim()[0],
             seq_matrix.get_dim()[1],
             motifs.get_dim()[0],
             n_iter,
             seq_matrix.get_dim()[1] - motifs.get_dim()[1] + 1,
             flip,
             n_threads),
      loglikelihood_max(n_row, 0.),
      cseq_layer(nullptr)
{
    this->loglikelihood_max = vector_d(n_row, 0.) ;

    // data and models
    // background motif (if any) is the last of the given motifs
    this->cseq_layer = new ConsensusSequenceLayer(std::move(seq_matrix),
                                                  std::move(motifs),
                                                  this->flip,
                                                  bckg_class) ;

    // intialise the class prob uniformly
    this->set_state_prob_uniform() ;
}


EMConsensusSequence::~EMConsensusSequence()
{   if(this->cseq_layer != nullptr)
    {   delete this->cseq_layer ;
        this->cseq_layer = nullptr ;
    }
    if(this->threads != nullptr)
    {   this->threads->join() ;
        delete this->threads ;
        this->threads = nullptr ;
    }
}

Matrix3D<double> EMConsensusSequence::get_sequence_models() const
{   return this->cseq_layer->get_model() ; }

EMConsensusSequence::exit_codes EMConsensusSequence::classify()
{
    size_t bar_update_n = this->n_iter ;
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
    return EMConsensusSequence::exit_codes::ITER_MAX ;
}

void EMConsensusSequence::compute_loglikelihood()
{   // compute the loglikelihood
    this->cseq_layer->compute_loglikelihoods(this->loglikelihood,
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
                                      std::bind(&EMConsensusSequence::compute_loglikelihood_routine,
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

void EMConsensusSequence::compute_loglikelihood_routine(size_t from,
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
                                     ConsensusSequenceLayer::p_min_log) ;
                }
            }
        }
    }
    done.set_value(true) ;
}

void EMConsensusSequence::compute_post_prob()
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
                                      std::bind(&EMConsensusSequence::compute_post_prob_routine,
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


void EMConsensusSequence::compute_post_prob_routine(size_t from,
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
                {
                    double p = std::max(this->post_prob(i,n_class,n_shift,n_flip) /
                                        this->post_prob_rowsum[i],
                                        ConsensusSequenceLayer::p_min) ;
                    this->post_prob(i,n_class,n_shift,n_flip) = p ;
                    colsums[n_class] += p ;
                }
            }
        }
    }
    post_prob_colsum.set_value(colsums) ;
}

void EMConsensusSequence::update_models()
{   this->cseq_layer->update_model(this->post_prob,
                                  this->threads) ;
}
