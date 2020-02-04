#ifndef CLASSSEQUENCEDATACREATOR_HPP
#define CLASSSEQUENCEDATACREATOR_HPP

#include <string>
#include <list>
#include <future>

#include <seqan/bam_io.h>  // BamFileIn
#include <seqan/bed_io.h>  // BedFileIn

#include <SequenceMatrixCreator.hpp>
#include <Matrix2D.hpp>
#include <Matrix3D.hpp>

/*!
 * \brief The ClassSequenceDataCreator class class allows to extract the data
 * that have been assigned to a given class, given a data partition.
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
class ClassSequenceDataCreator
{

    public:


        ClassSequenceDataCreator() = delete ;

        /*!
         * \brief Constructs an object to build a
         * class sequence matrix from a partition.
         * \param bed_file_path the path to the file containing
         * the references.
         * \param fasta_file_path the path to the file containing
         * the sequences.
         * \param prob_file_path the path to the file containing
         * the assignment probabilities of the partition.
         * It should be 4D matrix with the following dimensions :
         * 1st the number of regions, should be the number of
         * references in the BED file.
         * 2nd the number of classes.
         * 3rd the shifting freedom.
         * 4th the flipping freedom (1 for no flip, 2 otherwise).
         * \param from the upstream most relative position
         * to consider around the references. It may
         * be changed to make sure that the central bin
         * is centered on +/- 0.
         * \param to the dowmstream most relative position
         * to consider around the references. It may
         * be changed to make sure that the central bin
         * is centered on +/- 0.
         * \param class_k the index (1-based) of the class of
         * interest for which a matrix should be computed,
         * from the partition.
         */
        ClassSequenceDataCreator(const std::string& bed_file_path,
                                   const std::string& fasta_file_path,
                                   const std::string& prob_file_path,
                                   int from,
                                   int to,
                                   size_t class_k) ;

        /*!
         * Destructor.
         */
        ~ClassSequenceDataCreator() ;

        /*!
        * \brief Computes the matrix and returns it.
        * \return the class sequence matrix.
        * For each region, a consensus sequence is
        * returned as a probability matrix.
        * The matrix dimensions are :
        * 1st the number of regions.
        * 2nd the consensus sequence length.
        * 3rd 4 for A,C,G,T
        */
        Matrix3D<double> create_matrix() ;

    protected:
        /*!
         * \brief Bed file path.
         */
        std::string bed_file_path ;
        /*!
         * \brief Fasta file path.
         */
        std::string fasta_file_path ;
        /*!
         * \brief the path to the posterior probability
         * file (the partition).
         */
        std::string prob_file_path ;
        /*!
         * \brief The smallest relative coordinate from the region
         * center to consider (included).
         */
        int from ;
        /*!
         * \brief The biggest relative coordinate from the region
         * center to consider (not included).
         */
        int to ;
        /*!
         * \brief the class of interest (0-based).
         */
        size_t class_k ;

} ;

#endif // CLASSSEQUENCEDATACREATOR_HPP
