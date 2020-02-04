#include <EMBase.hpp>

#include <vector>
#include <stdexcept>   // std::invalid_argument
#include <future>      // std::promise, std::future
#include <utility>     // std::pair, std::move()
#include <functional>  // std::bind(), std::ref()
#include <numeric>     // std::iota()
#include <random>      // std::mt19937

#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>
#include <BetaDistribution.hpp>      // beta_distribution()
#include <Random.hpp>                // rand_string()
#include <RandomNumberGenerator.hpp> // getRandomNumberGenerator()
#include <Statistics.hpp>            // sd(), normal_pmf()


EMBase::EMBase(size_t n_row,
               size_t n_col,
               size_t n_class,
               size_t n_iter,
               size_t n_shift,
               bool flip,
               size_t n_threads=0)
    : n_row(n_row),
      n_col(n_col),
      n_class(n_class),
      n_shift(n_shift),
      flip(flip),
      n_flip(flip+1),
      n_iter(n_iter),
      l_model(n_col - n_shift + 1),
      loglikelihood(n_row, n_class, n_shift, n_flip, 0.),
      post_prob(n_row, n_class, n_shift, n_flip, 0.),
      post_state_prob(n_class, n_shift, n_flip, 0.),
      post_class_prob(n_class, 0.),
      post_prob_rowsum(n_row, 0.),
      post_prob_colsum(n_class, 0.),
      post_prob_tot(0.),
      threads(nullptr)
{
    // check n_shift value
    if(this->n_col < this->n_shift)
    {   char msg[4096] ;
        sprintf(msg, "Error! Shift is bigger than data column number "
                     "(%zu / %zu)!",
                this->n_shift, this->n_col) ;
        throw std::invalid_argument(msg) ;
    }
    /*
    // data structures
    this->loglikelihood = Matrix4D<double>(this->n_row,
                                           this->n_class,
                                           this->n_shift,
                                           this->n_flip,
                                           0.) ;
    this->post_prob = Matrix4D<double>(this->n_row,
                                       this->n_class,
                                       this->n_shift,
                                       this->n_flip,
                                       0.) ;
    this->post_state_prob = Matrix3D<double>(this->n_class,
                                             this->n_shift,
                                             this->n_flip,
                                             0.) ;
    this->post_class_prob = vector_d(this->n_class, 0) ;
    this->post_prob_rowsum = vector_d(this->n_row, 0) ;
    this->post_prob_colsum = vector_d(this->n_class, 0) ;
    this->post_prob_tot = 0 ;
    */
    // threads
    if(n_threads)
    {   this->threads = new ThreadPool(n_threads) ; }

}

EMBase::~EMBase()
{   // threads
    if(this->threads != nullptr)
    {   this->threads->join() ;
        delete this->threads ;
        this->threads = nullptr ;
    }
}

Matrix4D<double> EMBase::get_post_prob() const
{   return this->post_prob ; }

vector_d EMBase::get_post_class_prob() const
{   return this->post_class_prob ; }

void EMBase::set_state_prob_uniform()
{   double sum = this->n_class * this->n_shift * this->n_flip ;
    for(size_t i=0; i<this->n_class; i++)
    {   for(size_t j=0; j<this->n_shift; j++)
        {   for(size_t k=0; k<this->n_flip; k++)
            {   this->post_state_prob(i,j,k) = 1./sum ; }
        }
    }
}

void EMBase::set_post_prob_random(const std::string& seed)
{   // set random number generator
    // will be used to generate thread private seeds
    getRandomGenerator(seed) ;

    // don't parallelize
    if(this->threads == nullptr)
    {   std::promise<vector_d> promise ;
        std::future<vector_d> future = promise.get_future() ;
        this->set_post_prob_random_routine(0, this->n_row, seed, promise) ;
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
        // private seeds
        std::vector<std::string> private_seeds(n_threads) ;
        for(size_t i=0; i<n_threads; i++)
        {   futures[i] = promises[i].get_future() ;
            private_seeds[i] = rand_string(15) ;
        }

        // distribute work to threads
        // -------------------------- threads start --------------------------
        for(size_t i=0; i<n_threads; i++)
        {   // generate a private seed to set the random number generator
            // in this thread
            auto slice = slices[i] ;
            this->threads->addJob(std::move(
                                      std::bind(&EMBase::set_post_prob_random_routine,
                                                this,
                                                slice.first,
                                                slice.second,
                                                private_seeds[i],
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

    // compute class and state probs
    this->compute_class_prob() ;
}

void EMBase::set_post_prob_random_routine(size_t from,
                                          size_t to,
                                          const std::string& seed,
                                          std::promise<vector_d>& post_prob_colsum)
{   // random number generator
    std::mt19937 generator ;
    std::seed_seq seed_sequence(seed.begin(),seed.end()) ;
    generator.seed(seed_sequence) ;

    // this->post_prob_tot = 0. ;
    // this->post_prob_colsum = vector_d(this->n_class, 0.) ;
    vector_d colsums = vector_d(this->n_class, 0.) ;

    vector_d rowsums(this->n_row, 0) ;

    // random sampling
    beta_distribution<double> beta(1, this->n_row) ;
    for(size_t i=from; i<to; i++)
    {   for(size_t j=0; j<this->n_class; j++)
        {   for(size_t k=0; k<this->n_shift; k++)
            {   for(size_t l=0; l<this->n_flip; l++)
                {   double p = beta(generator) ;
                    this->post_prob(i,j,k,l) = p ;
                    rowsums[i] += p ;
                }
            }
        }
    }

    // normalization
    for(size_t i=from; i<to; i++)
    {   for(size_t j=0; j<this->n_class; j++)
        {   for(size_t k=0; k<this->n_shift; k++)
            {   for(size_t l=0; l<this->n_flip; l++)
                {   double p = this->post_prob(i,j,k,l) / rowsums[i] ;
                    this->post_prob(i,j,k,l) = p ;
                    // this->post_prob_tot        += p ;
                    // this->post_prob_colsum[j]  += p ;
                    colsums[j] += p ;
                }
            }
        }
    }

    // compute class and state probs
    // this->compute_class_prob() ;
    post_prob_colsum.set_value(colsums) ;
}

void EMBase::compute_class_prob()
{
    for(size_t n_class=0; n_class<this->n_class; n_class++)
    {   // reset total
        this->post_class_prob[n_class] = 0. ;
        for(size_t n_shift=0; n_shift<this->n_shift; n_shift++)
        {   for(size_t flip=0; flip<this->n_flip; flip++)
            {   // sum
                this->post_state_prob(n_class,n_shift,flip) = 0. ;
                for(size_t i=0; i<this->n_row; i++)
                {   this->post_state_prob(n_class,n_shift,flip) +=
                                                this->post_prob(i,n_class,n_shift,flip) ;
                }
                // normalize
                this->post_state_prob(n_class,n_shift,flip) /= this->post_prob_tot ;
                this->post_class_prob[n_class] += this->post_state_prob(n_class,n_shift,flip) ;
            }
        }
    }
}

void EMBase::center_post_state_prob()
{
    if(this->n_shift == 1)
    {   return ; }

    // the possible shift states
    vector_d shifts(this->n_shift) ;
    std::iota(shifts.begin(), shifts.end(), 1.) ;

    // the shift probabilities and the class probabilies
    // (no need to norm., class_prob sums to 1)
    double shifts_prob_measured_tot = 0. ;
    vector_d shifts_prob_measured(this->n_shift) ;
    for(size_t s=0; s<this->n_shift; s++)
    {   for(size_t k=0; k<this->n_class; k++)
        {   for(size_t f=0; f<this->n_flip; f++)
            {   shifts_prob_measured[s]  += this->post_state_prob(k,s,f) ;
                shifts_prob_measured_tot += this->post_state_prob(k,s,f) ;
            }
        }
    }


    // the shift mean and (biased) standard deviation
    double shifts_sd = sd(shifts, shifts_prob_measured, false) ;

    // the shift probabilities under the assumption that is
    // distributed as a gaussian centered on
    // the central shift state with sd and mean as in the data
    // sd as the data
    vector_d shifts_prob_centered(shifts.size(), 0.) ;
    double shifts_prob_centered_tot = 0. ;
    for(size_t i=0; i<shifts.size(); i++)
    {   shifts_prob_centered[i]   = normal_pmf(shifts[i],
                                               (this->n_shift/2)+1, shifts_sd) ;
        shifts_prob_centered_tot += shifts_prob_centered[i] ;
    }

    for(size_t k=0; k<this->n_class; k++)
    {   for(size_t f=0; f<this->n_flip; f++)
        {   for(size_t s=0; s<this->n_shift; s++)
            {   this->post_state_prob(k,s,f) = this->post_class_prob[k] *
                                                 shifts_prob_centered[s] /
                                                (this->n_flip * shifts_prob_centered_tot) ;
            }
        }
    }

    // shifts_prob_measured_tot = 0. ;
    shifts_prob_measured.clear() ;
    shifts_prob_measured.resize(this->n_shift) ;
    for(size_t s=0; s<this->n_shift; s++)
    {   for(size_t k=0; k<this->n_class; k++)
        {   for(size_t f=0; f<this->n_flip; f++)
            {   shifts_prob_measured[s]  +=
                        this->post_state_prob(k,s,f) ;
            }
        }
    }
}
