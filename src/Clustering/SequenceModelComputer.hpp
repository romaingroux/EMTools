#ifndef SEQUENCEMODELCOMPUTER_HPP
#define SEQUENCEMODELCOMPUTER_HPP

#include <ModelComputer.hpp>
#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

/*!
 * \brief The SequenceModelComputer class provides the
 * bases to compute the DNA sequence models, as computed
 * by one of the expectation maximization algorithm for
 * DNA sequence data, given a data matrix and a partitition.
 * The data partition must be a probability matrix, as returned
 * by one of the classesderiving EMBase.
 * The data must be a 2D matrix of dimensions NxL
 * containing the base at each of the Ldifferent positions in
 * N different regions.
 * The models returned are the weighted aggregations
 * of the sequences, assigned to each class.
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
class SequenceModelComputer : public ModelComputer
{
    public:

        /*!
        * \brief Constructs an object to retrieve
        * the sequence model given the data and their
        * classification results.
        * \param data the data.
        * Its dimensions are :
        * 1) the number of sequences
        * 2) the length of each sequence
        * The sequences should be encoded as follows :
        * A:0, C:1, G:2, T:3, else:5
        * \param post_prob the data class assignment
        * probabilities. Its dimensions are :
        * 1) the 1st dimension of the data matrix
        * 2) the number of classes
        * 3) the number of shift states
        * 4) the number of flip states
        * \param bckg_class  whether the last class of the
        * classification (posterior probabilities) is a
        * background class.
        * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
        */
        SequenceModelComputer(Matrix2D<int>&& data,
                              const Matrix4D<double>& post_prob,
                              bool bckg_class,
                              size_t n_threads) ;

        /*!
         * \brief Destructor.
         */
        virtual ~SequenceModelComputer() override ;

    protected:
        /*!
         * \brief the threads.
         */
        ThreadPool* threads ;

} ;

#endif // SEQUENCEMODELCOMPUTER_HPP
