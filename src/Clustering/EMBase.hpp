#ifndef EMBASE_HPP
#define EMBASE_HPP

#include <iostream>
#include <vector>
#include <future>  // std::promise

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>


typedef std::vector<double> vector_d ;


/*!
 * \brief The EMBase class is a base class
 * providing the basic support for classes
 * in implementing read density, sequence
 * and both at the time classification
 * procedures.
 */
class EMBase
{   public:
        /*!
         * \brief The possible exit codes for the classification
         * method.
         * 0 the classification procedure converged, 1 the
         * classification procedure ended by reaching the maximum
         * number of iterations, 2 the classification procedure
         * encountered an error.
         */
        enum exit_codes {CONVERGENCE=0, ITER_MAX, FAILURE} ;

    public:
         /*!
         * \brief Constructs an EMBase object.
         * \param n_row the number of rows in the data matrix.
         * \param n_col the number of columns in the data matrix.
         * \param n_iter the number of optimization iterations.
         * \param n_shift the number of shift states allowed.
         * \param flip whether flipping is allowed.
         * \param n_threads the number of parallel threads
         * to run the computations. 0 means no parallel
         * computing, everything is run on the main thread.
         * \throw std::invalid_argument if the shifting freedom
         * is bigger than the number of columns.
         */
        EMBase(size_t n_row,
               size_t n_col,
               size_t n_class,
               size_t n_iter,
               size_t n_shift,
               bool flip,
               size_t n_threads) ;

        EMBase(const EMBase& other) = delete ;

        /*!
         * \brief Destructor.
         */
        virtual ~EMBase() ;

        /*!
         * \brief Returns the posterior probability
         * of each point belonging to each class, for
         * each possible shift and flip state.
         * \return the posterior probability matrix,
         * with the following dimensions :
         * 1st dim : the number of region/sequence
         * 2nd dim : the number of classes
         * 3rd dim : the numbre of shift states
         * 4th dim : the number of flip states
         */
        virtual Matrix4D<double> get_post_prob() const ;

        /*!
         * \brief Returns the posterior class
         * probabilities (the total class
         * probability over all shift and
         * flip states).
         * \return the posterior class
         * probabilities.
         */
        virtual vector_d get_post_class_prob() const ;

        /*!
         * \brief Runs the models optimization and the
         * data classification.
         * \return a code indicating how the optimization
         * ended.
         */
        virtual EMBase::exit_codes classify() = 0 ;

    protected:

        /*!
         * \brief Computes the data log likelihood given the
         * current models, likelihood for each state.
         */
        virtual void compute_loglikelihood() = 0 ;

        /*!
         * \brief Computes the data posterior probabilties.
         */
        virtual void compute_post_prob() = 0 ;

        /*!
         * \brief Update the data models for all layers, given
         * the current posterior and class probabilities.
         */
        virtual void update_models() = 0;

        /*!
         * \brief Sets all the state probabilities
         * (all shift and flip states in all classes)
         * to a uniform probability.
         */
        void set_state_prob_uniform() ;

        /*!
         * \brief Sets the posterior
         * probabilities randomly (by
         * sampling them from a beta
         * distribution) and update all
         * other probabilities accordingly.
         * \param seed a seed to set the initial
         * state of the random number generator.
         */
        void set_post_prob_random(const std::string& seed) ;

        /*!
         * \brief The routine that effectively
         * sets the posterior probabilities randomly
         * (by sampling them from a beta
         * distribution).
         * \param from the index of the first row
         * in the data to consider.
         * \param to the index of the past last row
         * in the data to consider.
         * \param done the partial column (over the classes)
         * sum of posterior probabilities. If several routines
         * are running together, the colsums are retrieved by
         * summing up the vectors together.
         * Its length should be equal to the number of classes.
         * \param seed a seed to set the initial
         * state of the random number generator.
         */
        void set_post_prob_random_routine(size_t from,
                                          size_t to,
                                          const std::string& seed,
                                          std::promise<vector_d>& post_prob_colsum) ;

        /*!
         * \brief Computes the class/state probabilities from the
         * posterior probabilities.
         */
        void compute_class_prob() ;

        /*!
         * \brief Modifies the state probabilities in such a
         * way that the state probabilities are then normaly
         * distributed, centered on the middle shift state.
         * However, the overall class probabilities remain
         * unchanged.
         */
        void center_post_state_prob() ;

        /*!
         * \brief the number of rows in data.
         */
        size_t n_row ;
        /*!
         * \brief the number of columns in data.
         */
        size_t n_col ;
        /*!
         * \brief the number of classes.
         */
        size_t n_class ;
        /*!
         * \brief the number of shift states.
         */
        size_t n_shift ;
        /*!
         * \brief whther flip is allowed.
         */
        bool flip ;
        /*!
         * \brief zhe number of flip states.
         */
        size_t n_flip ;
        /*!
         * \brief the number of iterations.
         */
        size_t n_iter ;
        /*!
         * \brief the length of the models.
         */
        size_t l_model ;

        /*!
         * \brief the joint loglikelihood for each data point,
         * for each state (each class for each
         * shift and flip state). Its dimensions are :
         * 1st dim : the number of region/sequence
         * 2nd dim : the number of classes
         * 3rd dim : the numbre of shift states
         * 4th dim : the number of flip states
         */
        Matrix4D<double> loglikelihood ;
        /*!
         * \brief the posterior probabilities.
         * Its dimensions are :
         * 1st dim : the number of region/sequence
         * 2nd dim : the number of classes
         * 3rd dim : the numbre of shift states
         * 4th dim : the number of flip states
         * The values over a given region/sequence x
         * post_prob(x,.,.,.) should sum up to 1.0
         * and the entire matrix should sum up to
         * the 1st dimension value.
         */
        Matrix4D<double> post_prob ;
        /*!
         * \brief the states (shift and flip in each class)
         * probabilities.
         * Its dimensions are :
         * 1st dim : the number of classes
         * 2nd dim : the numbre of shift states
         * 3rd dim : the number of flip states
         * The values over a given class x
         * post_state_prob(x,.,.) should sum up to 1.0
         * and the entire matrix should sum up to
         * the 1st dimension value.
         */
        Matrix3D<double> post_state_prob ;
        /*!
         * \brief the total prob per class.
         * Its length should be equal to the
         * number of classes.
         */
        vector_d post_class_prob ;
        /*!
         * \brief the sum per row (data point) of post_prob.
         * Its length should be equal to the number of
         * regions/sequences.
         */
        vector_d post_prob_rowsum ;
        /*!
         * \brief the sum per column (class) of post_prob.
         * Its length should be equal to the number of
         * classes.
         */
        vector_d post_prob_colsum ;
        /*!
         * \brief the total of post_prob.
         */
        double post_prob_tot ;
        /*!
         * \brief the threads.
         */
        ThreadPool* threads ;
} ;


#endif // EMBASE_HPP
