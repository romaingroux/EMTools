#ifndef DATA2DLAYER_HPP
#define DATA2DLAYER_HPP

#include <DataLayer.hpp>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>

typedef std::vector<double> vector_d ;

/*!
 * \brief The Data2DLayer is a class that
 * implements the handling of probabilistic
 * models together with their data stored in
 * a 2D matrix.
 * A Data2DLayer is made of two parts :
 * 1) a 2D data matrix defined and implemented
 * here
 * 2) a model that is inherited from DataLayer
 * The models contains the paramters of a probabilistic
 * model with one or more classes that fits the data.
 * The data likelihood given the model can be computed
 * and the models can be uptaded given a set of posterior
 * probabilities representing the data assignment to the
 * different classes.
 */
class Data2DLayer : public DataLayer
{
    public:
        /*!
         * \brief Constructs an empty object.
         */
        Data2DLayer() ;

        /*!
         * \brief Constructs an object with the
         * given data.
         * An empty model is not initialised yet
         * as the model number of categories
         * depends on the final class.
         * \param data the data.
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        Data2DLayer(const Matrix2D<int>& data,
                    size_t n_class,
                    size_t n_shift,
                    bool flip,
                    bool bckg_class) ;

        /*!
         * \brief Constructs an object with the
         * given data.
         * An empty model is not initialised yet
         * as the model number of categories
         * depends on the final class.
         * \param data the data.
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        Data2DLayer(Matrix2D<int>&& data,
                    size_t n_class,
                    size_t n_shift,
                    bool flip,
                    bool bckg_class) ;

        /*!
         * \brief Constructs an object with the
         * given data and model.
         * The model dimensions set the number of
         * classes and the shifting freedom.
         * \param data the data.
         * \param the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
         Data2DLayer(const Matrix2D<int>& data,
                     const Matrix3D<double>& model,
                     bool flip,
                     bool bckg_class) ;


       /*!
        * \brief Constructs an object with the
        * given data and model.
        * The model dimensions set the number of
        * classes and the shifting freedom.
        * \param data the data.
        * \param the model.
        * \param flip whether flipping is allowed.
        * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
        */
        Data2DLayer(Matrix2D<int>&& data,
                    Matrix3D<double>&& model,
                    bool flip,
                    bool bckg_class) ;

        /*!
         * \brief Destructor.
         */
        virtual ~Data2DLayer() ;

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
        virtual void check_loglikelihood_dim(const Matrix4D<double>& loglikelihood) const override;

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
        virtual void check_loglikelihood_max_dim(const vector_d& loglikelihood_max) const override ;

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
        virtual void check_posterior_prob_dim(const Matrix4D<double>& posterior_prob) const override ;

        /*!
         * \brief the data.
         * Its dimensions are :
         * 1) the number of obervations / individuals
         * 2) the dimensionality of each individual
         * Example 1 : storing the number of reads mapped
         * along N different regions of lenght L would
         * require a matrix of dimensions NxL.
         * Example 2 : storing N different DNA sequences
         * of length L would require a matrix of dimensions
         * NxL. Because data is made to contains int, char
         * can be stored inside. However, for clarity, it may
         * be interesting to convert the DNA characters to
         * integer codes such as A:0, C:1, G:2, T:3
         */
        Matrix2D<int> data ;
        /*!
         * \brief the number of row in the data.
         */
        size_t n_row ;
        /*!
         * \brief the number of columns in the data.
         */
        size_t n_col ;

} ;

#endif // DATA2DLAYER_HPP
