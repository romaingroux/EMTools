#ifndef EMREADAPPLICATION_HPP
#define EMREADAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <string>

/*!
 * \brief The EMReadApplication class is a wrapper around an EMRead
 * instance creating an autonomous application to classify read/fragment
 * count data by directly  passing all the options and parameters from
 * the command line.
 */
class EMReadApplication: public ApplicationInterface
{
    public:
        EMReadApplication() = delete ;
        EMReadApplication(const EMReadApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        EMReadApplication(int argn, char** argv) ;

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
         * \brief the paths to the file containing the read
         * density data.
         */
        std::string file_read ;
        /*!
         * \brief the path to the file containing the
         * (0-based) index of the rows to filter out
         * from the data.
         */
        std::string file_filter ;
        /*!
         * \brief the path to the file in which the probability
         * matrix will be saved.
         */
        std::string file_out ;
        /*!
         * \brief the number of classes to partition the data into.
         */
        size_t n_class ;
        /*!
         * \brief the number of iterations allowed.
         */
        size_t n_iter ;
        /*!
         * \brief the shifting freedom.
         */
        size_t n_shift ;
        /*!
         * \brief whether flipping freedom is allowed.
         */
        bool flip ;
        /*!
         * \brief whether a constant class to model the
         * read count background should be added. This
         * class has mean number of read count at each
         * position.
         */
        bool bckg_class ;
        /*!
         * \brief the number of threads.
         */
        size_t n_threads ;
        /*!
         * \brief a seed to initialise the random number generator.
         */
        std::string seed ;
        /*!
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;
} ;

#endif // EMREADAPPLICATION_HPP
