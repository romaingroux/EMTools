#ifndef CLASSREADDATACREATOR_HPP
#define CLASSREADDATACREATOR_HPP

#include <string>
#include <list>
#include <future>

#include <seqan/bam_io.h>  // BamFileIn
#include <seqan/bed_io.h>  // BedFileIn

#include <CorrelationMatrixCreator.hpp>
#include <Matrix2D.hpp>

/*!
 * \brief The ClassReadDataCreator class allows to extract the data
 * that have been assigned to a given class, given a data partition.
 *
 * Given posterior probabilities and a read density matrix, the corresponding
 * class models can be computed. They are the weighted aggregations of the
 * data assigned to each given class. Instead of this, this program creates
 * the unfolded matrix that, if summed over the columns, gives the model of
 * class K.
 *
 * For a hard clustering methods, this procedure would simply correspond to the
 * creation of a matrix of dimensions N'xL where N'<=N is the number regions
 * assigned to class K among the N overall regions and L the number of bins of
 * the each region.
 *
 * In the case of a soft clustering methods, this procedure creates a matrix of
 * dimensions NxL' where L'=L-S+1 and S is the shifting freedom allowed during
 * the classification. The resulting matrix contains as many rows as the
 * starting matrix because in soft clustering, all regions (rows) are
 * assigned to all classes
 *
 * To construct a final matrix M3 of dimensions NxL3 where L3 covers a given
 * range <from>/<to>, the original matrix M1 of dimensions NxL is computed and
 * extended into a matrix M2 NxL2 with L2>=L1. The final M3 of dimensions NxL
 * is eventually computed, for class K, using the given posterior probabilities.
 * A row of the final matrix M3 is the weighted average of each of the S
 * possibles slices of the corresponding row in M2. The weights used are the
 * probabilities with which this row was assigned to class K, for each of
 * the S shift states, in each flip state.
 *
 * The original matrix M1 that was partitionned with shifting freedom S is
 * generated using the BED, BAM and BAI files that were originally used to
 * create it.
 * The posterior probabilities should be a 4D matrix in binary format, with
 * dimensions :
 * 1) number of regions
 * 2) number of classes
 * 3) number of shift states
 * 4) number of flip states
 * The results is returned through stdout as a 2D text matrix of dimensions :
 * 1) number of regions
 * 2) number of bins in each regions, as defined by the <from>/<to> range and
 * the <binsiz> value. ;
 */
class ClassReadDataCreator
{

    public:


        ClassReadDataCreator() = delete ;

        /*!
         * \brief Constructs an object to build a
         * class correlation matrix from a partition.
         * \param bed_file_path the path to the file containing
         * the references.
         * \param bam_file_path the path to the file containing
         * the targets.
         * \param bai_file_path the path to index file of the bam
         * file containing the targets.
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
         * \param bin_size the bin size in base pair.
         * \param class_k the index (1-based) of the class of
         * interest for which a matrix should be computed,
         * from the partition.
         * \param method how the targets should be counted.
         * READ all the positions inside the reads are
         * counted.
         * READ_ATAC only the +4bp position of +strand reads
         * and the -5bp of -strand reads are counted. It
         * correspond to the insertion position in ATAC-seq
         * data.
         * FRAGMENT all the positions within fragments (the
         * genome segment between a pair of reads, reads
         * included) are counted.
         * FRAGMENT_CENTER only the central position of the
         * fragements (the genome segment between a pair of
         * reads, reads included) are counted.
         */
        ClassReadDataCreator(const std::string& bed_file_path,
                                      const std::string& bam_file_path,
                                      const std::string& bai_file_path,
                                      const std::string& prob_file_path,
                                      int from,
                                      int to,
                                      int bin_size,
                                      size_t class_k,
                                      CorrelationMatrixCreator::methods method) ;

        /*!
         * Destructor.
         */
        ~ClassReadDataCreator() ;

        /*!
        * \brief Computes the matrix and returns it.
        * \return the class count matrix.
        * The weighted averages are floating values
        * however, they are rounded and returned as
        * integer values.
        */
        Matrix2D<int> create_matrix() ;


    protected:
        /*!
         * \brief Bed file path.
         */
        std::string bed_file_path ;
        /*!
         * \brief Bam file path.
         */
        std::string bam_file_path ;
        /*!
         * \brief Bam index file path.
         */
        std::string bai_file_path ;
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
         * \brief The bin size.
         */
        int bin_size ;
        /*!
         * \brief the class of interest (0-based).
         */
        size_t class_k ;
        /*!
         * \brief How to consider the sequenced fragments when computing
         * the bin values.
         */
        CorrelationMatrixCreator::methods method ;

} ;

#endif // CLASSREADDATACREATOR_HPP
