#ifndef SEQUENCELAYER_HPP
#define SEQUENCELAYER_HPP

#include <Data2DLayer.hpp>

#include <iostream>
#include <string>
#include <future>        // std::promise
#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

typedef std::vector<double> vector_d ;

/*!
 * \brief The SequenceLayer is a class that
 * implements the handling of probabilistic
 * models for DNA sequence data stored in
 * a 2D matrix. The characters should be
 * encoded using the following covention :
 * A:0, C:1, G:2, T:3, else:5-
 * A SequenceLayer is made of two parts :
 * 1) a 2D data matrix defined and implemented
 * in Data2DLayer
 * 2) a model that is inherited from Data2DLayer
 * but which behaviour is implemented here.
 * The data matrix dimensions are :
 * 1) the number of sequences
 * 2) the length of sequence
 * The data models are letter probability matrices
 * (1 per class) that contain the probability of
 * finding A, C, G or T at a given position in a
 * sequence. They are stored in a 3D matrix which
 * dimensions are :
 * 1) number of classes
 * 2) length of the model
 * 3) 4 for A, C, G, T respectively
 * model(0,x,0) + model(0,x,1) + model(0,x,2) +
 * model(0,x,3) must sum up to 1.0
 */
class SequenceLayer : public Data2DLayer
{
    public:
        /*!
         * \brief Computes the log-likelihood of the sub-
         * sequence - stored in a given row - and starting
         * at the offset <col> in the given sequence matrix.
         * The subsequence length is determined by the model
         * lenght.
         * \param seq the sequences in integer format.
         * Its dimensions should be :
         * 1) the number of sequences
         * 2) the length of the sequences
         * The following character encoding should be used :
         * A:0, C:1, G:2, T:3, else:5.
         * \param row the row containing the sequence of
         * interest.
         * \param col the index at which the sub-sequence
         * is starting (1st index inside the subsequence
         * of interest).
         * \param model_log a model containing the log
         * probability model. Its dimensions should be :
         * 1) the number of positions modeled
         * 2) 4 for A, C, G, T respectively
         * The columns should sum up to 1.
         * \return the log-likelihood of the sub-sequence
         * given the model.
         * \throw std::invalid_argument if 1) the offset is
         * invalid, 2) the sequence and the model have
         * incompatible dimensions or 3) the model 2n dimension
         * is not 4 (A,C,G,T).
         */
        static double score_subseq(const Matrix2D<int>& seq,
                                   size_t row,
                                   size_t col,
                                   const Matrix2D<double>& model_log) ;

    public:
        /*!
         * \brief Constructs an object with the
         * given data and an empty (0 values)
         * model.
         * \param data the data.
         * Its dimensions should be :
         * 1) the number of sequences
         * 2) the length of the sequences
         * Its dimensions should be :
         * 1) the number of sequences
         * 2) the length of the sequences
         * The following character encoding should be used :
         * A:0, C:1, G:2, T:3, else:5.
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        SequenceLayer(const Matrix2D<int>& data,
                      size_t n_class,
                      size_t n_shift,
                      bool flip,
                      bool bckg_class) ;

        /*!
         * \brief Constructs an object with the
         * given data and an empty (0 values)
         * model.
         * \param data the data.
         * Its dimensions should be :
         * 1) the number of sequences
         * 2) the length of the sequences
         * The following character encoding should be used :
         * A:0, C:1, G:2, T:3, else:5
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        SequenceLayer(Matrix2D<int>&& data,
                      size_t n_class,
                      size_t n_shift,
                      bool flip,
                      bool bckg_class) ;

        /*!
        * \brief Construct an object with the
        * given data and model.
        * The shifting freedom is set to (data number
        * of columns) - (the model 2nd dimension)
        * + 1.
        * \param data the data.
        * Its dimensions should be :
         * 1) the number of sequences
         * 2) the length of the sequences
         * The following character encoding should be used :
         * A:0, C:1, G:2, T:3, else:5
        * \param model the model with the following
        * dimensions :
        * 1) the number of classes
        * 2) the model length
        * 3) 4 for A, C, G, T respectively
        * model(x,y,0) + model(x,y,1) + model(x,y,2) +
        * model(x,y,3) should give 1.0
        * \param flip whether flipping is allowed.
        * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
        */
        SequenceLayer(const Matrix2D<int>& data,
                      const Matrix3D<double>& model,
                      bool flip,
                      bool bckg_class) ;


        /*!
        * \brief Construct an object with the
        * given data and model.
        * The shifting freedom is set to (data number
        * of columns) - (the model 2nd dimension)
        * + 1.
        * \param data the data.
        * Its dimensions should be :
         * 1) the number of sequences
         * 2) the length of the sequences
         * The following character encoding should be used :
         * A:0, C:1, G:2, T:3, else:5
        * \param model the model with the following
        * dimensions :
        * 1) the number of classes
        * 2) the model length
        * 3) 4 for A, C, G, T respectively
        * model(x,y,0) + model(x,y,1) + model(x,y,2) +
        * model(x,y,3) should give 1.0
        * \param flip whether flipping is allowed.
        * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
        */
        SequenceLayer(Matrix2D<int>&& data,
                      Matrix3D<double>&& model,
                      bool flip,
                      bool bckg_class) ;

        /*!
         * Destructor
         */
        virtual ~SequenceLayer() override ;

        /*!
         * \brief Computes the log likelihood of the data
         * given the current model parameters.
         * \param logliklihood a matrix to store the
         * results. It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param loglikelihood_max a vector containing the
         * max value for each row of loglikelihood.
         * Its length should be equal to the 1st dimension of
         * the loglikelihood matrix.
         * \throw std::invalid_argument if the dimensions are
         * incorrect.
         */
        virtual void compute_loglikelihoods(Matrix4D<double>& loglikelihood,
                                            vector_d& loglikelihood_max,
                                            ThreadPool* threads=nullptr) const override ;

        /*!
         * \brief Updates the model given the posterior
         * probabilities (the probabilities of each row
         * in the data to be assigned to each class,
         * for each shift and flip state).
         * It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * 0 values inside the model are forbidden.
         * The lowest possible value inside the model
         * is SequenceLayer::p_min.
         * \param posterior_prob the data assignment probabilities to
         * the different classes.
         */
        virtual void update_model(const Matrix4D<double>& posterior_prob,
                                  ThreadPool* threads=nullptr) override ;

    protected:
        /*!
         * \brief The routine that effectively performs the
         * loglikelihood computations.
         * \param from the index of the first row of the data
         * to considered.
         * \param to the index of the past last row of the data
         * to considered.
         * \param loglikelihood a matrix to store the
         * results. It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param loglikelihood_max a vector containing the
         * max value for each row of log_likelihood.
         * Its length should be equal to the 1st dimension of the
         * loglikelihood matrix.
         * \param model_log a vector containing the matrices with
         * the log values of the model for each class. Thus the
         * vector length should be equal to the 2nd dimension of
         * the loglikelihood matrix.
         * The matrix dimensions should be :
         * 1) the number of positions modeled
         * 2) 4 for A, C, G, T respectively
         * The columns should sum up to 1.
         * \param model_log_rev a vector containing the matrices with
         * the log values of the reverse strand model for each class
         * (the 1st position in the model becomes the last in the
         * reverse model and probabilities are swapped A<->T and C<->G).
         * Thus the vector length should be equal to the 2nd dimension
         * of the loglikelihood matrix.
         * The matrix dimensions should be :
         * 1) the number of positions modeled
         * 2) 4 for A, C, G, T respectively
         * The columns should sum up to 1.
         * \param done a promise to be filled when the routine
         * is done running.
         */
        void compute_loglikelihoods_routine(size_t from,
                                            size_t to,
                                            Matrix4D<double>& loglikelihood,
                                            vector_d& loglikelihood_max,
                                            const std::vector<Matrix2D<double>>& model_log,
                                            const std::vector<Matrix2D<double>>& model_log_rev,
                                            std::promise<bool>& done) const ;

        /*!
         * \brief The routine that effectively update the model.
         * 0 values inside the model are forbidden.
         * The lowest possible value inside the model
         * is SequenceLayer::p_min.
         * \param from the index of the first row of the
         * posterior probabilities to considered.
         * \param to the index of the past last row of the
         * posterior probabilities to considered.
         * \param posterior_prob the data assignment probabilities
         * to the different classes. It should have the following
         * dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param done a promise containing the partial model
         * computed from the given data slice. If several routines
         * work together at updating the model, they need to be
         * summed and normalized (by the column sum) to get the
         * final model.
         */
        void update_model_routine(size_t from,
                                  size_t to,
                                  const Matrix4D<double>& posterior_prob,
                                  std::promise<Matrix3D<double>>& done) const ;

        /*!
         * \brief Sets the last class as background class with the
         * mean base probababilities of the data at each position.
         */
        void create_bckg_class() ;

} ;
#endif // SEQUENCELAYER_HPP
