#ifndef READMODELCOMPUTER_HPP
#define READMODELCOMPUTER_HPP

#include <ModelComputer.hpp>
#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

/*!
 * \brief The ReadModelComputer class provides the
 * bases to compute the read density models, as computed
 * by one of the expectation maximization algorithm for
 * read density data, given a data matrix and a partitition.
 * The data partition must be a probability matrix, as
 * returned by one of the classes deriving EMBase.
 * The data must be a 2D matrix of dimensions NxL
 * containing the number of reads mapping at L
 * different positions in N different regions.
 * The models returned are the weighted aggregations
 * of the signal in each region, assigned to each
 * class.
 * The format of the model matrix returned is :
 * 1) the number of classes
 * 2) the model length + 1. The 1st column
 * contains the class probabilities. The signal
 * starts at the 2nd column. The 1st column sums
 * up to 1.0
 */
class ReadModelComputer : public ModelComputer
{
    public:

        /*!
        * \brief Constructs an object to retrieve
        * the read model given the data and their
        * classification results.
        * \param data the read density data.
        * Its dimensions are :
        * 1) the number of regions
        * 2) the length of each region
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
        ReadModelComputer(Matrix2D<int>&& data,
                          const Matrix4D<double>& post_prob,
                          bool bckg_class,
                          size_t n_threads) ;

        /*!
         * \brief Destructor.
         */
        virtual ~ReadModelComputer() override ;

    protected:
        /*!
         * \brief the threads.
         */
        ThreadPool* threads ;

} ;

#endif // READMODELCOMPUTER_HPP
