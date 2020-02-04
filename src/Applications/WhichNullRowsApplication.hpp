#ifndef WHICHNULLROWSAPPLICATION_HPP
#define WHICHNULLROWSAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <string>

/*!
 * \brief The WhichNullRowsApplication class implements an
 * autonomous application that reads a read density 2D matrix in
 * text format and returns the indices of the rows that contain
 * no signal (no read). The matrix dimensions should be :
 * 1) the number of regions
 * 2) the number of bins per region
 */
class WhichNullRowsApplication : public ApplicationInterface
{

    public:
        WhichNullRowsApplication()  = delete ;
        WhichNullRowsApplication(const WhichNullRowsApplication& other) = delete ;

        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        WhichNullRowsApplication(int argn, char** argv) ;

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
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;

} ;

#endif // WHICHNULLROWSAPPLICATION_HPP
