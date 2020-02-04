#ifndef CONSENSUSSEQUENCELAYER_HPP
#define CONSENSUSSEQUENCELAYER_HPP

#include <Data3DLayer.hpp>

#include <vector>
#include <future>   // std::promise

#include <ThreadPool.hpp>
#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>


typedef std::vector<double> vector_d ;

/*!
 * \brief The ConsensusSequenceLayer is a class that
 * implements the handling of probabilistic
 * models for DNA consensus sequence data stored in
 * a 3D matrix.
 * A DNA consensus matrix is a DNA sequence that
 * is represented using a letter probability
 * matrix (like the probabilistic model). Instead
 * of having A, C, G or T at each position, there
 * is aprobability distribution of each of these
 * nucleotides. Thus the data matrix contains a
 * collection of consensus sequences.
 * The data model is also a composed of letter
 * probability matrices.
 * The data matrix dimensions are :
 * 1) the number of sequences
 * 2) the length of sequence
 * 3) 4 for A, C, G and T respectively
 * A ConsusSequenceLayer is made of two parts :
 * 1) a 3D data matrix defined and implemented
 * in Data3DLayer
 * 2) a model that is inherited from Data3DLayer
 * but which behaviour is implemented here.
 */
class ConsensusSequenceLayer : public Data3DLayer
{
    public:
        /*!
         * \brief Computes the log-likelihood of the sub-
         * sequence - stored in a given row (1st dimension) -
         * and starting at the offset <col> (2nd dimension) in
         * the given sequence matrix.
         * The subsequence length is determined by the model
         * lenght.
         * \param cons_seq a probability matrix containing
         * the consensus sequence. Its dimensions are :
         * 1) the number of sequences
         * 2) the sequences length
         * 3) 4 for A,C,G,T.
         * cons_seq(x,y,0) + cons_seq(x,y,1) + cons_seq(x,y,2) +
         * cons_seq(x,y,3) should sum up to 1.0 as these
         * are the probabilities at a position y in sequence
         * x.
         * \param row the row containing the sequence of
         * interest.
         * \param col the index at which the sub-sequence
         * is starting (1st index inside the subsequence
         * of interest).
         * \param model a model containing the probability
         * model. Its dimensions are :
         * 1) the model length
         * 2) 4 for A, C, G, T respectively
         * The columns should sum up to 1.
         * \return the log-likelihood of the sub-sequence
         * given the model.
         * \throw std::invalid_argument if 1) the offset is
         * invalid, 2) the sequence and the model have
         * incompatible dimensions or 3) the model 2n dimension
         * is not 4 (A,C,G,T).
         */
        static double score_consensus_subseq(const Matrix3D<double>& cons_seq,
                                             size_t row,
                                             size_t col,
                                             const Matrix2D<double>& model) ;
    public:
        /*!
         * \brief Constructs an object with the
         * given data and an empty (0 values)
         * model.
         * \param data the data.
         * Its dimensions are :
         * 1) the number of sequences
         * 2) the sequences length
         * 3) 4 for A,C,G,T.
         * data(x,y,0) + data(x,y,1) + data(x,y,2) +
         * data(x,y,3) should sum up to 1.0 as these
         * are the probabilities at a position y in sequence
         * x.
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        ConsensusSequenceLayer(const Matrix3D<double>& data,
                               size_t n_class,
                               size_t n_shift,
                               bool flip,
                               bool bckg_class) ;

        /*!
         * \brief Constructs an object with the
         * given data and an empty (0 values)
         * model.
         * \param data the data.
         * Its dimensions are :
         * 1) the number of sequences
         * 2) the sequences length
         * 3) 4 for A,C,G,T.
         * data(x,y,0) + data(x,y,1) + data(x,y,2) +
         * data(x,y,3) should sum up to 1.0 as these
         * are the probabilities at a position y in sequence
         * x.
         * \param n_class the number of classes
         * of the model.
         * \param n_shift the number of shift
         * states of the model.
         * \param flip whether flipping is allowed.
         * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
         */
        ConsensusSequenceLayer(Matrix3D<double>&& data,
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
        * Its dimensions are :
        * 1) the number of sequences
        * 2) the sequences length
        * 3) 4 for A,C,G,T.
        * data(x,y,0) + data(x,y,1) + data(x,y,2) +
        * data(x,y,3) should sum up to 1.0 as these
        * are the probabilities at a position y in sequence
        * x.
        * \param model the model with the following
        * dimensions :
        * 1) the number of classes
        * 2) the model length
        * 3) 4 (A,C,G,T)
        * model(x,y,0) + model(x,y,1) + model(x,y,2) +
        * model(x,y,3) should sum up to 1.0 as this
        * are the probabilities of finding A, C, G, T
        * at position y in sequence x.
        * \param flip whether flipping is allowed.
        * \param bckg_class should the last class
        * of the model is constant and models the
        * background.
        */
        ConsensusSequenceLayer(const Matrix3D<double>& data,
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
        * Its dimensions are :
         * 1) the number of sequences
         * 2) the sequences length
         * 3) 4 for A,C,G,T.
         * data(x,y,0) + data(x,y,1) + data(x,y,2) +
         * data(x,y,3) should sum up to 1.0 as these
         * are the probabilities at a position y in sequence
         * x.
        * \param model the model with the following
        * dimensions :
        * 1) the number of classes
        * 2) the model length
        * 3) 4 (A,C,G,T)
        * model(x,y,0) + model(x,y,1) + model(x,y,2) +
        * model(x,y,3) should sum up to 1.0 as this
        * are the probabilities of finding A, C, G, T
        * at position y in sequence x.
        * \param flip whether flipping is allowed.
        * \param bckg_class should the last class
         * of the model is constant and models the
         * background.
        */
        ConsensusSequenceLayer(Matrix3D<double>&& data,
                               Matrix3D<double>&& model,
                               bool flip,
                               bool bckg_class) ;


        /*!
         * \brief Destructor.
         */
        virtual ~ConsensusSequenceLayer() ;

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
         * Its length should be equal to the loglikelihood matrix
         * 1st dimension.
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
         * is SequenceLayer::p_min.
         * \param posterior_prob the data assignment probabilities to
         * the different classes.
         * It should have the following dimensions :
         * 1st : same as the data number of row
         * 2nd : same as the model number of classes
         * 3rd : same as the number of shifts
         * 4th : same as the number of flip states
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
         * Its length should be equal to the loglikelihood 1st
         * dimension.
         * \param model_log a vector containing the matrices with
         * the log values of the model for each class.
         * The vector length should be equal to the number of classes
         * and each matrix dimensions should be :
         * 1) the model length
         * 2) 4 for A, C, G, T respectively
         * \param model_log_rev a vector containing the matrices with
         * the log values of the reverse strand model for each class
         * (the 1st position in the model becomes the last in the
         * reverse model and probabilities are swapped A<->T and C<->G).
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


    protected:




} ;


#endif // CONSENSUSSEQUENCELAYER_HPP
