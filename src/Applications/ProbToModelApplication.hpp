#ifndef PROBTOMODELAPPLICATION_HPP
#define PROBTOMODELAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <iostream>
#include <string>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>

/*!
 * \brief The ProbToModelApplication class is a wrapper around an ModelReferenceComputer
 * instance creating an autonomous application that computes data models given the
 * data and the results of the classification procedure (a 4D matrix stored in binary
 * posterior probability matrix of dimensions :
 * 1) number of regions
 * 2) number of classes
 * 3) number of shift states
 * 4) number of flip states)
 * The class models generated are stored in a 2D matrix in text format and
 * have the following dimensions :
 * 1) number of classes
 * 2) model length + 1, the 1st column contains the class probabilities
 */
class ProbToModelApplication : public ApplicationInterface
{
   public:
        ProbToModelApplication() = delete ;
        ProbToModelApplication(const ProbToModelApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        ProbToModelApplication(int argn, char** argv) ;
        /*!
         * \brief Destructor.
         */
        virtual ~ProbToModelApplication() override ;
        /*!
         * \brief Runs the application. The data model
         * is computed and displayed through the
         * stdout.
         * \return the exit code.
         */
        virtual int run() override ;

    private:
        /*!
         * \brief Parses the program command line options and
         * sets the object field accordingly.
         * If the help option is detected, the "runnable"
         * field is set to false and subsequent calls to
         * run() will produce nothing.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         * \throw std::invalid_argument if an error is found
         * in the program options.
         */
        void parseOptions(int argn, char** argv) ;

        /*!
         * \brief Checks that both matrices have
         * compatible dimensions.
         * For this, both matrices should have the
         * same number of rows (1st dimension) and
         * the data number of columns should not be
         * smaller than the probability matrix 3rd
         * dimension (shifting freedom incompatible
         * with data width).
         * \param data the data matrix.
         * \param prob the probability matrix.
         * \throw std::runtime_error with a descriptive
         * message if both matrices are not compatible.
         */
        void check_dimensions(const Matrix2D<int>& data,
                              const Matrix4D<double>& prob) const ;

        /*!
         * \brief Checks that both matrices have
         * compatible dimensions.
         * For this, both matrices should have the
         * 1st dimension size and the data 2nd
         * dimensions should not be smaller than
         * the probability matrix 3rd dimensions
         * (shifting freedom incompatible
         * with data width).
         * Additionaly, the data matrix should have
         * a 3rd dimension of size 4 (for A,C,G,T).
         * \param data the data matrix.
         * \param prob the probability matrix.
         * \throw std::runtime_error with a descriptive
         * message if both matrices are not compatible.
         */
        void check_dimensions(const Matrix3D<double>& data,
                              const Matrix4D<double>& prob) const ;

        /*!
         * \brief the path to the file containing the
         * read count data.
         */
        std::string file_read ;
        /*!
         * \brief the path to the file containing the
         * sequence data.
         */
        std::string file_seq ;
        /*!
         * \brief the path to the file containing the
         * consensus sequence data.
         */
        std::string file_consseq ;
        /*!
         * \brief the path to the file containing the
         * classification posterior probabilities.
         */
        std::string file_prob ;
        /*!
         * \brief the path to the file containing the
         * (0-based) index of the rows to filter out
         * from the data.
         */
        std::string file_filter ;
        /*!
         * \brief  whether the last class of the
         * classification (posterior probabilities) is a
         * background class.
         */
        bool bckg_class ;
        /*!
         * \brief the number of threads.
         */
        size_t n_threads ;
        /*!
         * \brief whether run() can be called.
         */
        bool runnable ;
} ;

#endif // PROBTOMODELAPPLICATION_HPP
