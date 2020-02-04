#ifndef MATRIXBINTOTXTAPPLICATION_HPP
#define MATRIXBINTOTXTAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <string>


/*!
 * \brief The MatrixBinToTxtApplication class implementes an autonomous
 * application that allows to convert a Matrix object dumped (using
 * Matrix.save(const std::string&) method) in binary format into text
 * format.
 * It can handle the conversion of Maitrx2D, 3D and 4D containing <int>,
 * <char> and <double> types.
 */
class MatrixBinToTxtApplication : public ApplicationInterface
{

    public:
        MatrixBinToTxtApplication()  = delete ;
        MatrixBinToTxtApplication(const MatrixBinToTxtApplication& other) = delete ;

        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        MatrixBinToTxtApplication(int argn, char** argv) ;

        /*!
         * \brief Runs the application. The data are classified
         * using the given settings and the posterior probability
         * matrix is returned through the stdout.
         * The matrix is a 4D matrix with dimensions :
         * regions, class, shift flip.
         * \return an exit code  EXIT_SUCCESS or EXIT_FAILURE
         * to return to the OS.
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
         * \brief the path to the file containing the matrix.
         */
        std::string file_matrix ;
        /*!
         * \brief the type of the data stored.
         */
        std::string type ;
        /*!
         * \brief the number of dimensions of the matrix.
         */
        size_t ndim ;
        /*!
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;

} ;

#endif // MATRIXBINTOTXTAPPLICATION_HPP
