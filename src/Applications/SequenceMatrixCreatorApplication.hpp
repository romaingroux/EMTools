#ifndef SEQUENCEMATRIXCREATORAPPLICATION_HPP
#define SEQUENCEMATRIXCREATORAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <string>

/*!
 * \brief The SequenceMatrixCreatorApplication class is a wrapper around a
 * SequenceMatrixCreator instance creating an autonomous application to
 * compute a DNA sequence matrix from a bed file containing the regions of
 * interest and a fasta file containing the genome sequences by directly
 * passing all the options and parameters from the command line.
 *
 * The sequence centers are defined as the center position of each
 * region contained in the bed file. The corresponding single bp
 * regions are then extended using the from/to parameters on each side.
 * The corresponding sequences are then extracted from the fasta file.
 *
 * The resulting matrix contains one sequence per row and one character per
 * column. DNA bases are converted into integer codes and printed on stdout.
 *
 */
class CorrelationMatrixCreatorApplication: public ApplicationInterface
{
    public:
        CorrelationMatrixCreatorApplication() = delete ;
        CorrelationMatrixCreatorApplication(const CorrelationMatrixCreatorApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        CorrelationMatrixCreatorApplication(int argn, char** argv) ;

        /*!
         * \brief Runs the application. A matrix containing the
         * sequences around the center of the bed regions, +/-
         * the given offsets is created by searching the fasta file
         * and printed on the stdout.
         * The matrix is a 2D matrix with dimensions :
         * 1) number of regions
         * 2) length of region (to - from + 1).
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
         * \brief the path to the bed file.
         */
        std::string file_bed ;
        /*!
         * \brief the path to the fasta file.
         */
        std::string file_fasta ;
        /*!
         * \brief a relative coordinate indicating the
         * most downstream position to consider around
         * each region in the bed file.
         */
        int from ;
        /*!
         * \brief a relative coordinate indicating the
         * most upstream position to consider around
         * each region in the bed file.
         */
        int to ;
        /*!
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;
} ;


#endif // SEQUENCEMATRIXCREATORAPPLICATION_HPP
