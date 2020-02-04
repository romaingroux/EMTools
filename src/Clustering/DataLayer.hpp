#ifndef DATALAYER_HPP
#define DATALAYER_HPP

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

typedef std::vector<double> vector_d ;

/*!
 * \brief The DataLayer is an abstract class that
 * defines the basic design to handle probabilistic
 * models together with their data.
 * A DataLayer is made of two parts :
 * 1) a data matrix (that should be defined and
 * implemented in children classes
 * 2) a model that is handled in this class
 * The model contains the parameters of a probabilistic
 * model with one or more classes that fits the data.
 * The data likelihood given the model can be computed
 * and the model can be updated given a set of
 * posterior probabilities representing the data
 * assignments to the different classes.
 */
class DataLayer
{
    public:
        /*!
         * \brief the smallest acceptable probability
         * for computations.
         */
        static const double p_min ;
        /*!
         * \brief the log of the smallest probability.
         */
        static const double p_min_log ;

        /*!
         * \brief The possible flip states.
         */
        enum flip_states{FORWARD=0, REVERSE} ;

        /*!
         * \brief Constructs an empty object.
         */
        DataLayer() ;

        /*!
         * \brief Constructs an object with the
         * given model parameters.
         * An empty model is not initialised yet
         * as the model number of categories
         * depends on the final class.
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        DataLayer(size_t n_class,
                  size_t n_shift,
                  bool flip,
                  bool bckg_class) ;

        /*!
         * \brief Constructs an object with the
         * given model.
         * The model dimensions set the number of
         * classes and the shifting freedom.
         * \param the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
         DataLayer(const Matrix3D<double>& model,
                   bool flip,
                   bool bckg_class) ;


       /*!
        * \brief Constructs an object with the
        * given model.
        * The model dimensions set the number of
        * classes and the shifting freedom.
        * \param the model.
        * \param flip whether flipping is allowed.
        * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
        */
        DataLayer(Matrix3D<double>&& model,
                  bool flip,
                  bool bckg_class) ;

        /*!
         * \brief Destructor.
         */
        virtual ~DataLayer() ;

        /*!
         * \brief Computes the log likelihood of the data
         * given the current model parameters.
         * \param loglikelihood a matrix to store the
         * results. It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param loglikelihood_max a vector containing the
         * max value for each row of log_likelihood.
         * Its length should be equal to the data row number.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        virtual void compute_loglikelihoods(Matrix4D<double>& loglikelihood,
                                            vector_d& loglikelihood_max,
                                            ThreadPool* threads=nullptr) const = 0 ;

        /*!
         * \brief Updates the model given the posterior
         * probabilities (the probabilities of each row
         * in the data to be assigned to each class,
         * for each shift and flip state).
         * \param posterior_prob the data assignment probabilities to
         * the different classes.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        virtual void update_model(const Matrix4D<double>& posterior_prob,
                                  ThreadPool* threads=nullptr) = 0 ;

        /*!
         * \brief Returns a copy of the current model.
         * \return the current model.
         * 1st dim : the number of classes
         * 2nd dim : the model length
         * 3rd dim : the number of value categories.
         */
        virtual Matrix3D<double> get_model() const ;

    protected:
        /*!
         * \brief Checks the argument has compatible
         * dimensions with the data and models. If this is
         * not the case, throw a std::invalid_argument with
         * a relevant message.
         * \param logliklihood a matrix to store the
         * results. It should have the following dimensions :
         * 1st : same as the data row number
         * 2nd : same as the model class number
         * 3rd : same as the shift state number
         * 4th : same as the flip state number
         * \throw std::invalid_argument if the dimensions are
         * incorrect.
         */
        virtual void check_loglikelihood_dim(const Matrix4D<double>& loglikelihood) const = 0;

        /*!
         * \brief Checks that the argument has compatible
         * dimensions with the data and models. If this is
         * not the case, throw a std::invalid_argument with
         * a relevant message.
         * \param loglikelihood_max a vector containing the
         * max value for each row of log_likelihood.
         * It should have a length equal to the number of
         * the data row number.
         * \throw std::invalid_argument if the dimensions are
         * incorrect.
         */
        virtual void check_loglikelihood_max_dim(const vector_d& loglikelihood_max) const = 0 ;

        /*!
         * \brief Checks the argument has compatible
         * dimensions with the data and models. If this is
         * not the case, throw a std::invalid_argument with
         * a relevant message.
         * \param posterior_prob a matrix to store the
         * results. It should have the following dimensions :
         * 1st : same as the data row number
         * 2nd : same as the model class number
         * 3rd : same as the shift state number
         * 4th : same as the flip state number
         * \throw std::invalid_argument if the dimensions are
         * incorrect.
         */
        virtual void check_posterior_prob_dim(const Matrix4D<double>& posterior_prob) const = 0 ;

        /*!
         * \brief the data model.
         * Its dimensions are :
         * 1) the number of classes
         * 2) the length of the model (the number of
         * positions/bins/values in a stretch of signal
         * that are represented by the model)
         * 3) the number of categories of values. For
         * continuous values, such as a number of reads
         * per position, this value is one. To handle
         * DNA sequences, this value would be 4 (A, C,
         * G, T) or even 5 to include 'N'.
         * Example 1 : representing 2 different classes
         * that models a read density signal of 100
         * bins of 1bp each mapping along a genome would
         * necessitate a matrix of dimensions 2x100x1.
         * Example 2 : representing 2 different classes of
         * DNA sequences of 12bp long as letter probability
         * matrices (A, C, G, T respectively) would require
         * a matrix of dimensions 2x12x4 and to ensure that
         * model(0,x,0) + model(0,x,1) + model(0,x,2) +
         * model(0,x,3) = 1.0 as it is a probability.
         */
        Matrix3D<double> model ;
        /*!
         * \brief whether flip is enabled.
         */
        bool flip ;
        /*!
         * \brief the number of classes in the model.
         */
        size_t n_class ;
        /*!
         * \brief the model length, its 2nd dimension.
         */
        size_t l_model ;
        /*!
         * \brief the number of variable categories in
         * the data. This is also the model 3rd
         * dimension.
         * Read counts are quantitative values and
         * have a number of categories equal to one
         * whereas as DNA sequences are made of
         * A,C,G,T (at least) and have 4 different
         * categories.
         */
        size_t n_category ;
        /*!
         * \brief the number of shift states.
         */
        size_t n_shift ;
        /*!
         * \brief the number of flip states.
         */
        size_t n_flip ;
        /*!
         * \brief A flag indicating that the last class of the model
         * is modelling the background. This class is considered constant
         * and should not be updated when calling update_model().
         */
        bool bckg_class ;
} ;

#endif // DATALAYER_HPP
