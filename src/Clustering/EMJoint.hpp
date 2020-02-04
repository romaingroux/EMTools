#ifndef EMJOINT_HPP
#define EMJOINT_HPP


#include <EMBase.hpp>

#include <vector>
#include <string>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ReadLayer.hpp>
#include <SequenceLayer.hpp>
#include <ConsensusSequenceLayer.hpp>


typedef std::vector<double> vector_d ;

/*!
 * \brief The EMJoint class provides a support to
 * handle the partition of regions based on one or
 * several read density/ies (called read layers)
 * together with (but optionally) the DNA sequences
 * of these regions.
 * This algorihm merges the partitioning algorithm
 * implemented by EMRead for the read density signal
 * handling and the partitioning algorithm implemented
 * by EMSequence for the DNA sequence partitioning.
 * Each class model thus contain one or more read
 * density model and (optionally) one DNA sequence
 * model (letter probability model).
 * For more informations about each algorithm, see
 * EMRead.hpp and EMSequence.hpp
 * For read density matrices, the dimensions are :
 * 1) the number of regions
 * 2) the region length
 * For the DNA sequence matrix, the dimensions are :
 * 1) the number of sequences
 * 2) the sequence lengths
 * The dimensions of all these matrices must be the
 * same. Additionally, the positions represented by
 * each cell of the read density and sequence matrices
 * should be striclty the same. For instance, if
 * read_matrix(1,1) contains the number of reads
 * mapping to position 123'434'123 from chromosome 1,
 * dna_matrix(1,1) should contain the base present
 * at the same postion. Thus, the bin size of the
 * read density matrices should be 1bp.
 * The probability matrix returned as results has
 * the following dimensions :
 * 1) the number of regions
 * 2) the number of classes
 * 3) the number of shift states
 * 4) the number of flip states
 * The sum over each row prob(x,.,.,.) should be
 * 1 and the sum over the entire matrix should be
 * equal to the 1st dimension value.
 */
class EMJoint : public EMBase
{
    public:
        /*!
         * \brief Constructs an object to partition the
         * region according to all the given read densities
         * with the given shifting and flipping freedom.
         * \param read_matrices a vector containing all
         * the different read density layers (ChIP-seq or related
         * signal) for the regions of interest. The matrix
         * dimensions are as follow :
         * 1) the number of regions,
         * 2) the region length.
         * \param n_class the number of classes to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model
         * the background by setting all its parameters, at all
         * positions, to the background base probabilties (for
         * sequence models) and the mean number of counts (for
         * read signal models). Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMJoint(const std::vector<Matrix2D<int>>& read_matrices,
                size_t n_class,
                size_t n_iter,
                size_t n_shift,
                bool flip,
                bool bckg_class,
                const std::string& seed="",
                size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * region according to all the given read densities
         * with the given shifting and flipping freedom.
         * \param read_matrices a vector containing all
         * the different read density layers (ChIP-seq or related
         * signal) for the regions of interest. The matrix
         * dimensions are as follow :
         * 1) the number of regions,
         * 2) the region length.
         * \param n_class the number of classes to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model
         * the background by setting all its parameters, at all
         * positions, to the background base probabilties (for
         * sequence models) and the mean number of counts (for
         * read signal models). Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMJoint(std::vector<Matrix2D<int>>&& read_matrices,
                size_t n_class,
                size_t n_iter,
                size_t n_shift,
                bool flip,
                bool bckg_class,
                const std::string& seed="",
                size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * region according to all the given read densities
         * and sequences with the given shifting and
         * flipping freedom.
         * \param read_matrices a vector containing all
         * the different read density layers (ChIP-seq or related
         * signal) for the regions of interest. The matrix
         * dimensions are as follow :
         * 1) the number of regions,
         * 2) the region length.
         * \param seq_matrix a matrix containing the DNA
         * sequences for the regions of interest. The matrix
         * dimensions are as follow:
         * 1) the number of sequences/regions,
         * 2) the sequence/region length.
         * \param n_class the number of classes to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model
         * the background by setting all its parameters, at all
         * positions, to the background base probabilties (for
         * sequence models) and the mean number of counts (for
         * read signal models). Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMJoint(const std::vector<Matrix2D<int>>& read_matrices,
                const Matrix2D<int>& seq_matrix,
                size_t n_class,
                size_t n_iter,
                size_t n_shift,
                bool flip,
                bool bckg_class,
                const std::string& seed="",
                size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * region according to all the given read densities
         * and sequences with the given shifting and
         * flipping freedom.
         * \param read_matrices a vector containing all
         * the different read density layers (ChIP-seq or related
         * signal) for the regions of interest. The matrix
         * dimensions are as follow :
         * 1) the number of regions,
         * 2) the region length.
         * \param seq_matrix a matrix containing the DNA
         * sequences for the regions of interest. The matrix
         * dimensions are as follow:
         * 1) the number of sequences/regions,
         * 2) the sequence/region length.
         * \param n_class the number of classes to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model
         * the background by setting all its parameters, at all
         * positions, to the background base probabilties (for
         * sequence models) and the mean number of counts (for
         * read signal models). Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMJoint(std::vector<Matrix2D<int>>&& read_matrices,
                Matrix2D<int>&& seq_matrix,
                size_t n_class,
                size_t n_iter,
                size_t n_shift,
                bool flip,
                bool bckg_class,
                const std::string& seed="",
                size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * region according to all the given read densities
         * consensus sequences with the given shifting and
         * flipping freedom.
         * \param read_matrices a vector containing all
         * the different read density layers (ChIP-seq or related
         * signal) for the regions of interest. The matrix
         * dimensions are as follow :
         * 1) the number of regions,
         * 2) the region length.
         * \param consseq_matrix a matrix containing the DNA
         * sequences represented as consensus sequences.
         * Instead of string sequences, the sequences are
         * represented as letter probability matrices.
         * Its dimensions are:
         * 1) the number of sequences/regions,
         * 2) the sequence/region length.
         * 3) 4 for A,C,G,T
         * consseq_matrix(x,y,0) + consseq_matrix(x,y,1) +
         * consseq_matrix(x,y,2) + consseq_matrix(x,y,3)
         * must sum up to 1.0 as these are the probabilities
         * of A, C, G, T at position y of sequence x.
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model
         * the background by setting all its parameters, at all
         * positions, to the background base probabilties (for
         * sequence models) and the mean number of counts (for
         * read signal models). Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMJoint(const std::vector<Matrix2D<int>>& read_matrices,
                const Matrix3D<double>& consseq_matrix,
                size_t n_class,
                size_t n_iter,
                size_t n_shift,
                bool flip,
                bool bckg_class,
                const std::string& seed="",
                size_t n_threads=0) ;

        /*!
         * \brief Constructs an object to partition the
         * region according to all the given read densities
         * consensus sequences with the given shifting and
         * flipping freedom.
         * \param read_matrices a vector containing all
         * the different read density layers (ChIP-seq or related
         * signal) for the regions of interest. The matrix
         * dimensions are as follow :
         * 1) the number of regions,
         * 2) the region length.
         * \param consseq_matrix a matrix containing the DNA
         * sequences represented as consensus sequences.
         * Instead of string sequences, the sequences are
         * represented as letter probability matrices.
         * Its dimensions are:
         * 1) the number of sequences/regions,
         * 2) the sequence/region length.
         * 3) 4 for A,C,G,T
         * consseq_matrix(x,y,0) + consseq_matrix(x,y,1) +
         * consseq_matrix(x,y,2) + consseq_matrix(x,y,3)
         * must sum up to 1.0 as these are the probabilities
         * of A, C, G, T at position y of sequence x.
         * \param n_class the number of region classes
         * to search.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param bckg_class the last class is used to model
         * the background by setting all its parameters, at all
         * positions, to the background base probabilties (for
         * sequence models) and the mean number of counts (for
         * read signal models). Since the background is constant,
         * this class will never be updated.
         * \param seed a seed to initialise the random number
         * generator.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         */
        EMJoint(std::vector<Matrix2D<int>>&& read_matrices,
                Matrix3D<double>&& consseq_matrix,
                size_t n_class,
                size_t n_iter,
                size_t n_shift,
                bool flip,
                bool bckg_class,
                const std::string& seed="",
                size_t n_threads=0) ;

        EMJoint(const EMJoint& other) = delete ;

        /*!
         * \brief Destructor.
         */
        virtual ~EMJoint() override ;

        /*!
         * \brief Returns all layer read models.
         * The models are in the same order
         * as the data were given to the
         * constructor.
         * \return a vector containing the
         * models.
         * The matrix dimensions are :
         * 1) the number of classes
         * 2) the length of the model
         * 3) 1
         * The vector lenght is equal to the number
         * of read layers.
         */
        std::vector<Matrix3D<double>> get_read_models() const ;

        /*!
         * \brief Returns the sequence models.
         * \return a vector containing the
         * models.
         * The matrix has the following dimensions :
         * 1) the number of classes
         * 2) the model length
         * 3) 4 for A,C,G,T
         * motif(x,y,0) + motif(x,y,1) + motif(x,y,2) +
         * motif(x,y,3) sum up to one as these
         * corresponds to the probabilies of A, C, G, T
         * at position y of model x.
         */
        Matrix3D<double> get_sequence_models() const ;

        /*!
         * \brief Returns the consensus sequence
         * models.
         * \return a vector containing the
         * models.
         * The matrix has the following dimensions :
         * 1) the number of classes
         * 2) the model length
         * 3) 4 for A,C,G,T
         * motif(x,y,0) + motif(x,y,1) + motif(x,y,2) +
         * motif(x,y,3) sum up to one as these
         * corresponds to the probabilies of A, C, G, T
         * at position y of model x.
         */
        Matrix3D<double> get_consensus_sequence_models() const ;

        /*!
         * \brief Runs the sequence model optimization and
         * the data classification.
         * \return a code indicating how the optimization
         * ended.
         */
        virtual EMJoint::exit_codes classify() override ;

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
         * \brief This is a routine of compute_loglikelihood() that
         * computes the joint loglikelihood by summing the
         * individual loglikelihood obtained from each data layer.
         * At the same time, this method rescales the loglikelihood
         * values by substacting to each value the maximum
         * loglikelihood value found in the same data row,
         * for each layer.
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
         * It length should be the number of classes.
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
         * \brief the number of data layers.
         */
        size_t n_layer ;
        /*!
         * \brief the log likelihood buffers for each individual
         * layer (one element per layer).
         * Each matrix dimensions are :
         * 1) dim : the number of region/sequence
         * 2) dim : the number of classes
         * 3) dim : the numbre of shift states
         * 4) dim : the number of flip states
         */
        std::vector<Matrix4D<double>> loglikelihood_layer ;
        /*!
         * \brief the max loglikelihood value for
         * each each data layer (1st dimension)
         * and each data row of the given layer
         * (2nd dimension).
         * Its length should be equal to the
         * number of data layers and each item
         * (vector inside the vector) length
         * should be equal to the number of
         * regions in the data.
         */
        std::vector<vector_d> loglikelihood_max ;

        /*!
         * \brief A vector containing the pointers
         * to the objects managing all the read
         * data and models.
         */
        std::vector<ReadLayer*> read_layers ;
        /*!
         * \brief A pointer to the object managing
         * the sequence data and their models.
         */
        SequenceLayer* seq_layer ;
        /*!
         * \brief A pointer to the object managing
         * the consensus sequences data and their
         * models.
         */
        ConsensusSequenceLayer* consseq_layer ;
} ;

#endif // EMJOINT_HPP
