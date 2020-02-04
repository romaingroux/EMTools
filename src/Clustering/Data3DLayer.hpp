#ifndef DATA3DLAYER_HPP
#define DATA3DLAYER_HPP

#include <DataLayer.hpp>

#include <Matrix3D.hpp>
#include <Matrix4D.hpp>

typedef std::vector<double> vector_d ;

class Data3DLayer : public DataLayer
{
    public:
        /*!
         * \brief Constructs an empty object.
         */
        Data3DLayer() ;

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
        Data3DLayer(const Matrix3D<double>& data,
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
        Data3DLayer(Matrix3D<double>&& data,
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
         Data3DLayer(const Matrix3D<double>& data,
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
        Data3DLayer(Matrix3D<double>&& data,
                    Matrix3D<double>&& model,
                    bool flip,
                    bool bckg_class) ;

        /*!
         * \brief Destructor.
         */
        virtual ~Data3DLayer() ;


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
         * 1) the number of observations / individuals
         * 2) the dimensionality of each individual
         * 3) the number of categories possible for each
         * dimension of a given indivual.
         * Example 1 : storing N DNA sequences of length L
         * as probability matrices, instead of a regular
         * string sequences, would require a matrix of
         * dimensions NxLx4 (A, C, G, T respectively). The
         * 1st position of the 1st sequence would be
         * stored in
         * probability of A : data(0,0,0)
         * probability of C : data(0,0,1)
         * probability of G : data(0,0,2)
         * probability of T : data(0,0,3)
         * It should be ensured that these values sum up
         * to 1 because these are probabilities.
         */
        Matrix3D<double> data ;
        /*!
         * \brief The size of data 1st dimension.
         */
        size_t n_dim1 ;
        /*!
         * \brief The size of data 2nd dimension.
         */
        size_t n_dim2 ;
        /*!
         * \brief The size of data 3rd dimension.
         */
        size_t n_dim3 ;

} ;

#endif // DATA3DLAYER_HPP
