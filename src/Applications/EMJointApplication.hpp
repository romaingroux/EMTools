#ifndef EMJOINTAPPLICATION_HPP
#define EMJOINTAPPLICATION_HPP

#include <ApplicationInterface.hpp>
#include <EMJoint.hpp>

#include <string>

/*!
 * \brief The EMJointApplication class is a wrapper around an EMJoint
 * instance creating an autonomous application to classify one or more
 * read/fragment cout data and/or DNA sequences by directly
 * passing all the options and parameters from the command line.
 */
class EMJointApplication: public ApplicationInterface
{
    public:
        EMJointApplication() = delete ;
        EMJointApplication(const EMJointApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        EMJointApplication(int argn, char** argv) ;

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
         * \brief a coma separated list of paths to the files
         * containing the read density data
         */
        std::string files_read ;
        /*!
         * \brief the path to the file containing the
         * sequence data.
         */
        std::string file_sequence ;
        /*!
         * \brief the path to the file containing the
         * consensus sequence data.
         */
        std::string file_conssequence ;
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
         * sequence background should be added. This
         * class has the sequence background probabilities
         * (for the sequence model) and the mean number of
         * read counts (for read signal model) at each
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


#endif // EMJOINTAPPLICATION_HPP
