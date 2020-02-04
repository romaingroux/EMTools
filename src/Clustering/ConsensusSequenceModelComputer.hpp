#ifndef CONSENSUSSEQUENCEMODELCOMPUTER_HPP
#define CONSENSUSSEQUENCEMODELCOMPUTER_HPP

#include <ModelComputer.hpp>
#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

/*!
 * \brief The ConsensusSequenceModelComputer class provides the
 * bases to compute the DNA sequence models, as computed
 * by one of the expectation maximization algorithm for
 * consensus DNA sequence data, given a data matrix and a partitition.
 * The data partition must be a probability matrix, as returned
 * by one of the classes deriving EMBase.
 * The data must be a 3D matrix of dimensions NxLx4
 * containing the probability of each base at each of the L different
 * positions in N different regions.
 * The models returned are the weighted aggregations
 * of the consensus sequences, assigned to each class.
 * The format of the model matrix returned is :
 * 1) 4 * the number of classes
 * 2) the model length + 1.
 * The rows should by read 4 by 4 to get each class
 * model. For instance, if the partition contains 3
 * classes, the model matrix has 12 rows. The 1st
 * class is contained in rows 1-4, the 2nd class in
 * rows 5-8 and the 3rd class in rows 9-12.
 * The 1st column always contains the class probabilities. The signal
 * starts at the 2nd column. For each class, the 1st column sums
 * up to 1.0
 */
class ConsensusSequenceModelComputer : public ModelComputer
{
    public:

        /*!
        * \brief Constructs an object to retrieve
        * the (consensus) sequence model given the data
        * and their classification results.
        * \param data the data.
        * Its dimensions are :
        * 1) the number of consensus sequences
        * 2) the length of each consensus sequence
        * 3) 4 for A, C, G, T
        * \param post_prob the data class assignment
        * probabilities. Its dimensions are :
        * 1) the 1st dimension of the data matrix
        * 2) the number of classes
        * 3) the number of shift states
        * 4) the number of flip states
        * \param bckg_class whether the last class of the
        * classification (posterior probabilities) is a
        * background class.
        * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
        */
        ConsensusSequenceModelComputer(Matrix3D<double>&& data,
                              const Matrix4D<double>& post_prob,
                              bool bckg_class,
                              size_t n_threads) ;

        /*!
         * \brief Destructor.
         */
        virtual ~ConsensusSequenceModelComputer() override ;

    protected:
        /*!
         * \brief the threads.
         */
        ThreadPool* threads ;

} ;

#endif // CONSENSUSSEQUENCEMODELCOMPUTER_HPP
