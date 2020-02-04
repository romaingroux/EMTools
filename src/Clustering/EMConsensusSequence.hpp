#ifndef EMCONSENSUSSEQUENCE_HPP
#define EMCONSENSUSSEQUENCE_HPP

#include <EMBase.hpp>

#include <vector>
#include <string>
#include <future>       // std::promise

#include <Matrix3D.hpp>
#include <ConsensusSequenceLayer.hpp>


typedef std::vector<double> vector_d ;

/*!
 * \brief The EMConsensusSequence class provides a support to handle
 * read density data and to cluster them using a probabilistic
 * clustering algorithm. This algorithm is a modification of
 * ChIPPartitioning to cluster DNA consensus sequences instead of
 * read density signal (for ChIPPartitioning see
 * https://academic.oup.com/bioinformatics/article/30/17/2406/2748187).
 * DNA consensus sequences are a representation of DNA sequences
 * using letter probability matrices. At each position, instead of
 * having either A, C, G or T, there is a probability distribution
 * of these 4 nucleotides.
 *
 * In brief, the consensus sequences are stored in a NxLx4 dimension
 * matrix that contains the probability of the bases at each of the
 * L positions of N sequences.
 * The base probabilities at each position of each consensus
 * sequence are modeled as having being sampled from K different
 * A,C,G,T distributions (where K is the number of
 * classes).
 * The K classes are represented as K letter probability
 * matrices of dimensions L-S+1x4 that represent the
 * probability of sampling A, C, G, T at a given position
 * and where S is the number of shift states.
 * To allow the algorithm to search DNA motifs (class
 * models) shorter than the entire consensus sequences, S shift
 * states are allowed.
 * This results in the creation of S overlapping slices of
 * length L-S+1 for each sequence. If S=1, then no
 * shifting is allowed and the slices correspond to
 * the entire consensus sequences (L-1+1 = L).
 * Each slice is then compared with the K different
 * models, which have the same length, and a probability
 * (a similarity score) is computed.
 * To better the data realignment process, flipping
 * the slices can also be performed (which correspond
 * to searching the DNA motif on the reverse strand).
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
 * 1) the number of sequences
 * 2) the number of classes
 * 3) the number of shift states
 * 4) the number of flip states
 * The sum over each row prob(x,.,.,.) should be
 * 1 and the sum over the entire matrix should be
 * equal to the 1st dimension value.
 */
class EMConsensusSequence : public EMBase
{
    public:
        /*!
         * \brief Constructs an object to partition the
         * given consensus sequences (rows) according to
         * their motif content.
         * The sequences models are initialised randomly.
         * \param sequence_matrix a matrix containing the
         * consensus sequences in a probability matrix.
         * Its dimensions are :
         * 1) the number of consensus sequences
         * 2) the length of the consensus sequences
         * 3) 4 for A,C,G,T
         * The sums over the 1st and 2nd dimensions should
         * be 1. The overall sum of the matrix values should
         * be the st dimension.
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model the background
         * by setting all its parameters, at all positions, to the
         * background base probabilties. Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMConsensusSequence(const Matrix3D<double>& sequence_matrix,
                            size_t n_class,
                            size_t n_iter,
                            size_t n_shift,
                            bool flip,
                            bool bckg_class,
                            const std::string& seed="",
                            size_t n_threads=0) ;
        /*!
         * \brief Constructs an object to partition the
         * given consensus sequences (rows) according to
         * their motif
         * content.
         * The sequences models are initialised randomly.
         * \param sequence_matrix a matrix containing the
         * consensus sequences in a probability matrix.
         * Its dimensions are :
         * 1st the number of consensus sequences
         * 2nd the length of the consensus sequences
         * 3rd 4 for A,C,G,T
         * The sums over the 1st and 2nd dimensions should
         * be 1. The overall sum of the matrix values should
         * be the st dimension.
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model the background
         * by setting all its parameters, at all positions, to the
         * background base probabilties. Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMConsensusSequence(Matrix3D<double>&& sequence_matrix,
                            size_t n_class,
                            size_t n_iter,
                            size_t n_shift,
                            bool flip,
                            bool bckg_class,
                            const std::string& seed="",
                            size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * given consensus sequences (rows) according to
         * their motif content.
         * The sequences class models are initialised using
         * the given motifs. The class probabilities are
         * initialised uniformelly.
         * The shifting freedom is set to (data  2n dimension)
         * - (the model 2nd dimension) + 1.
         * \param a matrix containing the consensus sequences
         * in a probability matrix. Its dimensions are :
         * 1st the number of consensus sequences
         * 2nd the length of the consensus sequences
         * 3rd 4 for A,C,G,T
         * The sums over the 1st and 2nd dimensions should
         * be 1. The overall sum of the matrix values should
         * be the st dimension.
         * \param motifs a matrix containing the different initial
         * class models with the following dimensions :
         * dim1 the number of classes
         * dim2 the model length
         * dim3 4 for A,C,G,T
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param flip whether flipping is allowed.
         * \param bckg_class indicates that the last class in the
         * given motifs is used to model the background and it
         * should never be updated.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMConsensusSequence(const Matrix3D<double>& sequence_matrix,
                            const Matrix3D<double>& motifs,
                            size_t n_iter,
                            bool flip,
                            bool bckg_class,
                            size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * given consensus sequences (rows) according to
         * their motif content.
         * The sequences class models are initialised using
         * the given motifs. The class probabilities are
         * initialised uniformelly.
         * The shifting freedom is set to (data  2n dimension)
         * - (the model 2nd dimension) + 1.
         * \param a matrix containing the consensus sequences
         * in a probability matrix. Its dimensions are :
         * 1st the number of consensus sequences
         * 2nd the length of the consensus sequences
         * 3rd 4 for A,C,G,T
         * The sums over the 1st and 2nd dimensions should
         * be 1. The overall sum of the matrix values should
         * be the st dimension.
         * \param motifs a matrix containing the different initial
         * class models with the following dimensions :
         * dim1 the number of classes
         * dim2 the model length
         * dim3 4 for A,C,G,T
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param flip whether flipping is allowed.
         * \param bckg_class indicates that the last class in the
         * given motifs is used to model the background and it
         * should never be updated.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMConsensusSequence(Matrix3D<double>&& sequence_matrix,
                            Matrix3D<double>&& motifs,
                            size_t n_iter,
                            bool flip,
                            bool bckg_class,
                            size_t n_threads=0) ;

        EMConsensusSequence(const EMConsensusSequence& other) = delete ;

        /*!
         * \brief Destructor.
         */
        virtual ~EMConsensusSequence() override ;

        /*!
         * \brief Returns the class sequence model.
         * \return the class sequence model.
         */
        Matrix3D<double> get_sequence_models() const ;

        /*!
         * \brief Runs the sequence model optimization and
         * the data classification.
         * \return a code indicating how the optimization
         * ended.
         */
        virtual EMConsensusSequence::exit_codes classify() override ;

    private:

        /*!
         * \brief Computes the data log likelihood given the
         * current models, for all layers and the joint
         * likelihood for each state (the sum of the layer
         * likelihoods for all layers, for a given state).
         * To avoid numerical issues when computing posterior
         * probabilities, the lowest possible value authorized
         * as log likelihood is ConsensusSequenceLayer::p_min_log.
         * Any value below is replaced by this one.
         */
        virtual void compute_loglikelihood() override ;

        /*!
         * \brief This is a routine of compute_loglikelihood().
         * This method rescales the loglikelihood values by
         * substacting to each value the maximum loglikelihood
         * value found in the same data row.
         * To avoid numerical issues when computing posterior
         * probabilities, the lowest possible value authorized
         * as log likelihood is ConsensusSequenceLayer::p_min_log.
         * Any value below is replaced by this one.
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
         * ConsensusSequenceLayer::p_min. Any value below
         * is replaced by this one.
         */
        virtual void compute_post_prob() override ;

        /*!
         * \brief The routine that effectively computes
         * the posterior probabilties.
         * To avoid numerical issues the lowest possible
         * value authorized as posterior probability is
         * ConsensusSequenceLayer::p_min. Any value below
         * is replaced by this one.
         * \param from the index of the first row
         * in the data to consider.
         * \param to the index of the past last row
         * in the data to consider.
         * \param done the partial column (over the classes)
         * sum of posterior probabilities. If several routines
         * are running together, the colsums are retrieved by
         * summing up the vectors together.
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
        ConsensusSequenceLayer* cseq_layer ;

} ;

#endif // EMCONSENSUSSEQUENCE_HPP
