#ifndef READLAYER_HPP
#define READLAYER_HPP

#include <Data2DLayer.hpp>

#include <future>       // std::promise

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>
#include <ThreadPool.hpp>

typedef std::vector<double> vector_d ;

/*!
 * \brief The ReadLayer is a class that
 * implements the handling of probabilistic
 * models for read density data stored in
 * a 2D matrix.
 * A ReadLayer is made of two parts :
 * 1) a 2D data matrix defined and implemented
 * in Data2DLayer
 * 2) a model that is inherited from Data2DLayer
 * but which behaviour is implemented here.
 * The data matrix dimensions are :
 * 1) the number of regions
 * 2) the length of each region
 * The data model is stored in a 3D matrix
 * which dimensions are :
 * 1) number of classes
 * 2) length of the models
 * 3) 1
 * For a given class, each value is the lambda
 * parameter of a Poisson distribution that
 * models the distribution of the number of
 * reads that are observed in regions, at
 * the corresponding position.
 */
class ReadLayer : public Data2DLayer
{
    public:
        /*!
         * \brief Constructs an object with the
         * given data and an empty (0 values)
         * model.
         * \param data the data.
         * Its dimensions should be :
         * 1) the number of regions
         * 2) the length of the regions
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        ReadLayer(const Matrix2D<int>& data,
                  size_t n_class,
                  size_t n_shift,
                  bool flip,
                  bool bckg_class,
                  ThreadPool* threads = nullptr) ;

        /*!
         * \brief Constructs an object with the
         * given data and an empty (0 values)
         * model.
         * \param data the data.
         * Its dimensions should be :
         * 1) the number of regions
         * 2) the length of the regions
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        ReadLayer(Matrix2D<int>&& data,
                  size_t n_class,
                  size_t n_shift,
                  bool flip,
                  bool bckg_class,
                  ThreadPool* threads = nullptr) ;

        /*!
         * \brief Construct an object with the
         * given data and model.
         * \param data the data.
         * Its dimensions should be :
         * 1) the number of regions
         * 2) the length of the regions
         * \param the model.
         * Its dimensions are :
         * 1) the number of classes
         * 2) the number of positions modeled
         * 3) 1
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        ReadLayer(const Matrix2D<int>& data,
                  const Matrix3D<double>& model,
                  bool flip,
                  bool bckg_class,
                  ThreadPool* threads = nullptr) ;

        /*!
         * \brief Construct an object with the
         * given data and model.
         * \param data the data.
         * Its dimensions should be :
         * 1) the number of regions
         * 2) the length of the regions
         * \param the model.
         * Its dimensions are :
         * 1) the number of classes
         * 2) the number of positions modeled
         * 3) 1
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        ReadLayer(Matrix2D<int>&& data,
                  Matrix3D<double>&& model,
                  bool flip,
                  bool bckg_class,
                  ThreadPool* threads = nullptr) ;

        /*!
         * Destructor
         */
        virtual ~ReadLayer() override ;

        /*!
         * \brief Computes the log likelihood of the data
         * given the current model parameters.
         * During this process, a normalized version of the
         * models, having a sum of signal of 1 count in average,
         * is used (a copy of the models is normalized, meaning
         * that the original models can still be retrieved the
         * dedicated getter).
         * \param logliklihood a matrix to store the
         * results. It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param loglikelihood_max a vector containing the
         * max value for each row of loglikelihood.
         * Its length should be equal to the data row number.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
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
         * 0 values inside the model are forbidden.
         * The lowest possible value inside the model
         * is ReadLayer::p_min.
         * \param posterior_prob the data assignment
         * probabilities to the different classes.
         * It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        virtual void update_model(const Matrix4D<double>& posterior_prob,
                                  ThreadPool* threads=nullptr) override ;

        /*!
         * \brief Updates the model given the posterior
         * probabilities (the probabilities of each row
         * in the data to be assigned to each class,
         * for each shift and flip state).
         * This method does the same as the virtual method it
         * overloads. The only difference is that, for run time
         * gain, it is given the sum over the columns of the
         * posterior_prob matrix which is computed by the virtual
         * method.
         * 0 values inside the model are forbidden.
         * The lowest possible value inside the model
         * is ReadLayer::p_min.
         * \param posterior_prob the data assignment probabilities to
         * the different classes.
         * It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
         * \param posterior_prob_colsum the sum over the columns
         * (classes) of the posterior_prob matrix. Thus its length
         * should be equal to the 2nd dimension of the matrix.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        void update_model(const Matrix4D<double>& posterior_prob,
                          const vector_d& posterior_prob_colsum,
                          ThreadPool* threads=nullptr) ;
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
         *
         * \param done a promise to be filled Its length should be equal to the 1st dimension of the
         * loglikelihood matrix.when the routine
         * is done running.
         */
        void compute_loglikelihoods_routine(size_t from,
                                                    size_t to,
                                                    Matrix4D<double>& loglikelihood,
                                                    vector_d& loglikelihood_max,
                                                    std::promise<bool>& done) const ;

        /*!
         * \brief The routine that effectively update the model.
         *  0 values inside the model are forbidden.
         * The lowest possible value inside the model
         * is ReadLayer::p_min.
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
         * \param posterior_prob_colsum the sum over the columns
         * (classes) of the posterior_prob matrix. Thus its length
         * should be equal to the 2nd dimension of the matrix.
         * \param promise a promise containing the partial model
         * computed from the given data slice. If several routines
         * work together to update the model, the promise matrices
         * need to be summed up to get the final model.
         */
        void update_model_routine(size_t from,
                                  size_t to,
                                  const Matrix4D<double>& posterior_prob,
                                  const vector_d& posterior_prob_colsum,
                                  std::promise<Matrix3D<double>>& promise) const ;

        /*!
         * \brief Computes the mean number of reads present in
         * each slice (of length l_model), in each row
         * of the data and store them in this->window_means.
         * \param threads a pointer to a thread pool to
         * parallelize the computations. If nullptr is given,
         * the computations are performed by the main thread.
         */
        void compute_window_means(ThreadPool* threads) ;

        /*!
         * \brief The routine that effectively computes the
         * window means.
         * \param from the index of the first row of the
         * data to considered.
         * \param to the index of the past last row of the
         * data to considered.
         * \param done a promise to fill when the routine
         * is done running.
         */
        void compute_window_means_routine(size_t from,
                                          size_t to,
                                          std::promise<bool>& done) ;

        /*!
         * \brief Checks that the argument has compatible
         * dimensions with the data and models. If this is
         * not the case, throw a std::invalid_argument with
         * a relevant message.
         * \param posterior_class_prob a vector containing the
         * class probabilities.
         * It should have a length equal to the number of
         * classes.
         * \throw std::invalid_argument if the dimensions are
         * incorrect.
         */
        void check_posterior_prob_colsum_dim(const vector_d& posterior_prob_colsum) const ;

        /*!
         * \brief Sets the last class as background class with the mean
         * number of counts in the data at each position.
         */
        void create_bckg_class() ;

        /*!
         * \brief contains the data means, for
         * each window of size l_model.
         */
        Matrix2D<double> window_means ;

} ;

#endif // READLAYER_HPP
