#ifndef CORRELATIONMATRIXCREATORAPPLICATION_HPP
#define CORRELATIONMATRIXCREATORAPPLICATION_HPP

#include <ApplicationInterface.hpp>
#include <CorrelationMatrixCreator.hpp>  // CorrelationMatrixCreator::methods

#include <string>

/*!
 * \brief The CorrelationMatrixCreatorApplication class is a wrapper around a
 * CorrelationMatrixCreator instance creating an autonomous application to
 * compute a read/fragment count matrix from a BAM file by directly passing
 * all the options and parameters from the command line.
 *
 * The center of each region is computed as the center of the BED regions.
 * Then a set of equally sized, non-overlapping bins, centered on the region
 * center and covering the interval [from,to] is build for each region.
 * Then, each bin is assigned the number of read/fragment positions (targets)
 * present in the BAM file that are mapped at that position.
 *
 * The read/fragment counts are then computed, for each bin,
 * from the BAM file.
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
         * read densities around the center of the bed regions, +/-
         * the given offsets is created by searching the BAM file
         * and printed on the stdout.
         * The matrix is a 2D matrix with dimensions :
         * 1) number of regions
         * 2) length of region (to - from + 1) / bin_size.
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
         * \brief the path to the bam file.
         */
        std::string file_bam ;
        /*!
         * \brief the path to the bam index file.
         */
        std::string file_bai ;
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
         * \brief the size of the bin that will be used
         * to bin the signal in the regions [from,to] around
         * each region in the bed file.
         */
        int bin_size ;
        /*!
         * \brief How to consider the sequenced fragments when computing
         * the bin values.
         */
        CorrelationMatrixCreator::methods method ;
        /*!
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;
} ;


#endif // CORRELATIONMATRIXCREATORAPPLICATION_HPP
