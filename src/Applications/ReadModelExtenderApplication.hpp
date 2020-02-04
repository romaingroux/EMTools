#ifndef READMODELEXTENDERAPPLICATION_HPP
#define READMODELEXTENDERAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <iostream>
#include <string>

#include <CorrelationMatrixCreator.hpp>

/*!
 * \brief The ReadModelExtenderApplication class is a class implementing an
 * autonomous application that extends read models of length L' (L' = L - S + 1
 * where L is the number of column of the data matrix and S the
 * shifting freedom allowed during the classification) to a new model
 * length L'' = L' + E  (where E is the number of columns to add to the
 * models) given the data matrix and the results of the classification
 * (posterior probability matrix).
 * To do this, the read density matrix from which the original models
 * were computed is generated and extended by adding E/2 columns on each side.
 * Eventually the extended models are computed from the extended matrix using 
 * the given posterior probabities.
 * The posterior probabilies should be stored in a 4D matrix in binary format
 * of dimensions :
 * 1) number of regions
 * 2) number of classes
 * 3) number of shift states
 * 4) number of flip states
 * The extended model are returned through the stdout as a 2D text matrix of
 * dimensions :
 * 1) number of classes
 * 2) L''+1, the 1st column contains the class probabilities.
 */
class ReadModelExtenderApplication : public ApplicationInterface
{
   public:
        ReadModelExtenderApplication() = delete ;
        ReadModelExtenderApplication(const ReadModelExtenderApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        ReadModelExtenderApplication(int argn, char** argv) ;
        /*!
         * \brief Destructor.
         */
        virtual ~ReadModelExtenderApplication() override ;
        /*!
         * \brief Runs the application. The data new model
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
         * \brief the path to the file containing the
         * classification posterior probabilities.
         */
        std::string file_prob ;
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
         * \brief the number of columns to add to the
         * matrix (half of this value on each side).
         */
        int ext ;
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

#endif // READMODELEXTENDERAPPLICATION_HPP
