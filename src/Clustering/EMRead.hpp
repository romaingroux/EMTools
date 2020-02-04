#ifndef EMREAD_HPP
#define EMREAD_HPP

#include <EMBase.hpp>

#include <vector>
#include <string>
#include <future>       // std::promise

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <ReadLayer.hpp>


typedef std::vector<double> vector_d ;

/*!
 * \brief The EMRead class provides a support to handle
 * read density data and to cluster them using a probabilistic
 * clustering algorithm (called ChIPPartitioning, see
 * https://academic.oup.com/bioinformatics/article/30/17/2406/2748187)
 * In brief, the data are stored in a NxL dimension matrix
 * that contains the number of reads that map at
 * L different positions along N regions.
 * The signal at each position of each region is
 * modelled as having being sampled from K different
 * Poisson distributions (where K is the number of
 * classes).
 * The K classes are represented as K vectors of L-S+1
 * values that are the lambda parameters of the Poisson
 * distributions and S is the number of shift states.
 * To allow the algorithm to realign a region and a
 * class model, S shift states are allowed. This
 * result in the creation of S overlapping slices of
 * length L-S+1 for each region. If S=1, then no
 * shifting is allowed and the slices correspond to
 * the entire regions (L-1+1 = L).
 * Each slice is then compared with the K different
 * models, which have the same length, and a probability
 * (a similarity score) is computed.
 * To better the data realignment process, flipping
 * the slices can also be performed.
 * Once all slices have been compared with all models,
 * and that the probabilities have been computed, the
 * class models are updated using a weighted ungapped
 * alignment. Each class model is the weighted
 * aggregation of all data slices. The weights are the
 * probability with which each slice as been
 * assigned to that class.
 * This scheme is repeated iteratively to perform
 * an expectation-maximization optimization of the
 * partition.
 * This probability matrix dimensions are :
 * 1) the number of regions
 * 2) the number of classes
 * 3) the number of shift states
 * 4) the number of flip states
 * The sum over each row prob(x,.,.,.) should be
 * 1 and the sum over the entire matrix should be
 * equal to the 1st dimension value.
 */
class EMRead : public EMBase
{
    public:
        /*!
         * \brief Constructs an object to partition the
         * regions (rows) according to the shape of the signal
         * with the given shifting and flipping freedom.
         * \param read_matrix a matrix containing the read
         * densitiy (ChIP-seq or related signal) for the
         * regions of interest.
         * Its dimensions are :
         * 1) the number of regions
         * 2) the length of the regions
         * \param n_class the number of region classes
         * to search.Matrix
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model the
         * background by setting all its parameters, at all
         * positions, to the mean number of counts. Since
         * the background is constant, this class will never
         * be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMRead(const Matrix2D<int>& read_matrix,
               size_t n_class,
               size_t n_iter,
               size_t n_shift,
               bool flip,
               bool bckg_class,
               const std::string& seed="",
               size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * region (rows) according to the shape of the signal
         * with the given shifting and flipping freedom.
         * \param read_matrix a matrix containing the read
         * densitiy (ChIP-seq or related signal) for the
         * regions of interest.
         * Its dimensions are :
         * 1) the number of regions
         * 2) the length of the regions
         * \param n_class the number of region classes
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model the
         * background by setting all its parameters, at all
         * positions, to the mean number of counts. Since
         * the background is constant, this class will never
         * be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMRead(Matrix2D<int>&& read_matrix,
               size_t n_class,
               size_t n_iter,
               size_t n_shift,
               bool flip,
               bool bckg_class,
               const std::string& seed="",
               size_t n_threads=0) ;

        EMRead(const EMRead& other) = delete ;

        /*!
         * \brief Destructor.
         */
        virtual ~EMRead() override ;

        /*!
         * \brief Returns the class read signal model.
         * \return the class read signal model.
         * Its dimensions are :
         * 1) the number of classes
         * 2) the length of the model
         * 3) 1
         * \param n_class the number of region classes
         */
        Matrix3D<double> get_read_models() const ;

        /*!
         * \brief Runs the read signal model optimization and
         * the data classification.
         * \return a code indicating how the optimization
         * ended.
         */
        virtual EMRead::exit_codes classify() override ;

    private:

        /*!
         * \brief Computes the data log likelihood given the
         * current models, for all layers and the joint
         * likelihood for each state (the sum of the layer
         * likelihoods for all layers, for a given state).
         * To avoid numerical issues when computing posterior
         * probabilities, the lowest possible value authorized
         * as log likelihood is ReadLayer::p_min_log. Any
         * value below is replaced by this one.
         */
        virtual void compute_loglikelihood() override ;

        /*!
         * \brief This is a routine of compute_loglikelihood().
         * This method rescales the loglikelihood values by
         * substacting to each value the maximum loglikelihood
         * value found in the same data row.
         * To avoid numerical issues when computing posterior
         * probabilities, the lowest possible value authorized
         * as log likelihood is ReadLayer::p_min_log. Any
         * value below is replaced by this one.
         * \param from the index of the first row
         * in the data to consider.
         * \param to the index of the past last row
         * in the data to consider.
         * \param done a promise to fill when the method
         * is done.
         */
        void compute_loglikelihood_routine(size_t from,
                                           size_t to,
                                           std::promise<bool>& done) ;

        /*!
         * \brief Computes the data posterior probabilties.
         * To avoid numerical issues the lowest possible
         * value authorized as posterior probability is
         * ReadLayer::p_min. Any value below is replaced
         * by this one.
         */
        virtual void compute_post_prob() override ;

        /*!
         * \brief The routine that effectively computes
         * the posterior probabilties.
         * To avoid numerical issues the lowest possible
         * value authorized as posterior probability is
         * ReadLayer::p_min. Any value below is replaced
         * by this one.
         * \param from the index of the first row
         * in the data to consider.
         * \param to the index of the past last row
         * in the data to consider.
         * \param done the partial column (over the classes)
         * sum of posterior probabilities. If several routines
         * are running together, the colsums are retrieved by
         * summing up the vectors together.
         * Its length should be the number of classes.
         */
        void compute_post_prob_routine(size_t from,
                                       size_t to,
                                       std::promise<vector_d>& post_prob_colsum) ;

        /*!
         * \brief Update the data models for all layers, given
         * the current posterior and class probabilities.
         */
        virtual void update_models() override ;

        /*!
         * \brief the max loglikelihood value for
         * each data row.
         */
        std::vector<double> loglikelihood_max ;
        /*!
         * \brief A pointer to the object managing
         * the data and their model.
         */
        ReadLayer* read_layer ;

} ;

#endif // EMREAD_HPP
