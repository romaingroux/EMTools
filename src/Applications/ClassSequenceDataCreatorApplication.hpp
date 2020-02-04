#ifndef CLASSSEQUENCEDATACREATORAPPLICATION_HPP
#define CLASSSEQUENCEDATACREATORAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <string>

#include <SequenceMatrixCreator.hpp>  // SequenceMatrixCreator::methods

/*!
 * \brief The ClassSequenceDataCreatorApplication is a wrapper
 * around a ClassSequenceDataCreator that creates a sequence
 * matrix containing the data assigned to a given class, given
 * a partition of these data.
 *
 * Given posterior probabilities and a sequence matrix, the corresponding
 * class models can be computed. They are the weighted aggregations of the
 * DNA sequences assigned to each given class. However, because DNA sequences
 * cannot be summed, the aggregation are represented as probability matrices
 * or consensus sequence (A+C is represented as 50%A, 50%C, 0%G, 0%T). Instead
 * of this, this program creates the unfolded matrix that, if summed over the
 * columns, gives the model of class K.
 *
 * For a hard clustering methods, this procedure would simply correspond to the
 * creation of a matrix of dimensions N'xL where N'<=N is the number sequences
 * assigned to class K among the N overall sequences and L the length of
 * the each sequence.
 *
 * In the case of a soft clustering methods, this procedure creates a 3D matrix of
 * dimensions NxL'x4. This matrix contains N probability matrices, each one of
 * dimensions L'x4 where L'=L-S+1, 4 corresponds to A, C, G, T and S is the
 * shifting freedom allowed during the classification. The resulting matrix
 * contains as many rows as the starting matrix because in soft clustering, all
 * sequences (rows) are assigned to all classes
 *
 * To construct a final matrix M3 of dimensions NxL3 where L3 covers a given
 * range <from>/<to>, the original matrix M1 of dimensions NxL is computed and
 * extended into a matrix M2 NxL2 with L2>=L1. The final M3 of dimensions NxL
 * is eventually computed, for class K, using the given posterior probabilities.
 * A row of the final matrix M3 is the weighted average of each of the S
 * possibles slices of the corresponding row in M2, represented as a probability
 * matrix. The weights used are the probabilities with which this row was assigned
 * to class K, for each of the S shift states, in each flip state.
 *
 * The original matrix M1 that was partitionned with shifting freedom S is
 * generated using the BED and fasta files that were originally used to
 * create it.
 * The posterior probabilities should be a 4D matrix in binary format, with
 * dimensions :
 * 1) number of sequences
 * 2) number of classes
 * 3) number of shift states
 * 4) number of flip states
 * The results is returned as a 3D binary matrix of dimensions :
 * 1) number of sequences
 * 2) length of the sequences, as defined by the <from>/<to> range
 * 3= 4 for A, C, G, T
 */
class ClassSequenceDataCreatorApplication : public ApplicationInterface
{
    public:
        ClassSequenceDataCreatorApplication() = delete ;
        ClassSequenceDataCreatorApplication(
                const ClassSequenceDataCreatorApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        ClassSequenceDataCreatorApplication(int argn, char** argv) ;

        /*!
         * \brief TODO
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
         * \brief the path to the posterior probability
         * file (the partition).
         */
        std::string file_prob ;
        /*!
         * \brief the path to the file in which the
         * results will be written.
         */
        std::string file_out ;
        /*!
         * \brief the coordinate of the most upstream
         * position that was in the original matrix, in
         * relative coordinate.
         */
        int from ;
        /*!
         * \brief the coordinate of the most downstream
         * position that was in the original matrix, in
         * relative coordinate.
         */
        int to ;
        /*!
         * \brief the class of interest (0-based).
         */
        size_t class_k ;
        /*!
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;
} ;

#endif // CLASSSEQUENCEDATACREATORAPPLICATION_HPP
