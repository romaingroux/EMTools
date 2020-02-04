#ifndef CORRELATIONMATRIXCREATOR_HPP
#define CORRELATIONMATRIXCREATOR_HPP

#include <string>
#include <list>
#include <future>

#include <seqan/bam_io.h>  // BamFileIn
#include <seqan/bed_io.h>  // BedFileIn

#include <ReadMatrixCreator.hpp>
#include <Matrix2D.hpp>

/*!
 * \brief The CorrelationMatrixCreator class allows
 * to create read density correlation matrices.
 * A read density correlation matrix contains the number
 * of reads/fragments mapped at different positions around
 * a set of reference positions.
 * This class reads the reference positions from
 * a BED file and the targets from a BAM file, with the
 * help of an BAI file to speed up the parsing. For each
 * reference, the region center is computed and then a
 * region covering the interval [from,to] is build
 * around the middle and divided into equally sized
 * bins. Finally, each bin is assigned the number of
 * targets present in the BAM file that are mapped at
 * that position.
 * The final matrix contains one row per reference,
 * with the  number of targets counted at each possible
 * position (bin). relative to this reference.
 */
class CorrelationMatrixCreator: public ReadMatrixCreator
{
    public:

        CorrelationMatrixCreator() = delete ;

        /*!
         * \brief Constructs an object to build a
         * correlation matrix.
         * \param bed_file_path the path to the file containing
         * the references.
         * \param bam_file_path the path to the file containing
         * the targets.
         * \param bai_file_path the path to index file of the bam
         * file containing the targets.
         * \param from the upstream most relative position
         * to consider around the references. It may
         * be changed to make sure that the central bin
         * is centered on +/- 0.
         * \param to the dowmstream most relative position
         * to consider around the references. It may
         * be changed to make sure that the central bin
         * is centered on +/- 0.
         * \param bin_size the bin size in base pair.
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
        CorrelationMatrixCreator(const std::string& bed_file_path,
                                 const std::string& bam_file_path,
                                 const std::string& bai_file_path,
                                 int from,
                                 int to,
                                 int bin_size,
                                 CorrelationMatrixCreator::methods method) ;
        /*!
         * Destructor.
         */
        virtual ~CorrelationMatrixCreator() ;

        /*!
        * \brief Computes the matrix and returns it.
        * \return the count matrix.
        */
        virtual Matrix2D<int> create_matrix() override ;

    protected:
        /*!
         * \brief Seek in the BAM file right before the last
         * record upstream the given region. The margin
         * parameters allows to modify the region start
         * value.
         * To read a record within the region, a read
         * operation is required to get ride of the
         * record right
         * \param region the region in front of which the
         * pointer is desired.
         * \param margin
         * which streams in the stream vectors to use.
         * \return whether the reading pointer could be moved
         * to the desired position.
         */
        bool jump_upstream(const GenomeRegion& region,
                           int margin) ;

        /*!
         * \brief A generic routine that reads the following records
         * until finding the first one located downstream the region
         * of interest (the definition of the first target downstream
         * the region of interest depends if READ/READ_ATAC/FRAGMENT
         * or FRAGMENT_CENTER is set as method).
         * All record overlapping the region of interest are stored
         * in the target lists.
         * The reading pointer is supposed to be located
         * upstream the region of interest. If this is note the case,
         * the method will read records until reaching the end of
         * the file.
         * \param region the region of interest.
         */
        void to_downstream_target(const GenomeRegion& region) ;

        /*!
         * \brief The routine that reads the following records
         * until finding the first one located downstream the region
         * of interest if READ or READ_ATAC is set as method.
         * All record overlapping the region of interest are stored
         * in the target lists.
         * The reading pointer is supposed to be located
         * upstream the region of interest. If this is note the case,
         * the method will read records until reaching the end of
         * the file.
         * \param region the region of interest.
         */
        void to_downstream_read(const GenomeRegion& region) ;

        /*!
         * \brief The routine that reads the following records
         * until finding the first one located downstream the region
         * of interest if FRAGMENT or FRAGMENT_CENTER is set as
         * method.
         * All record overlapping the region of interest are stored
         * in the target lists.
         * The reading pointer is supposed to be located
         * upstream the region of interest. If this is note the case,
         * the method will read records until reaching the end of
         * the file.
         * \param region the region of interest.
         */
        void to_downstream_fragment(const GenomeRegion& region) ;

        /*!
         * \brief Clear the content of the target lists.
         */
        void clear_target_lists() ;

        /*!
         * \brief Update the given row of the count matrix with
         * the content of the target lists.
         * \param matrix_row_index the index of the row, in the
         * count matrix.
         */
        void update_count_matrix(size_t row_index) ;

        /*!
         * \brief A buffers containing the
         * target mapped on the forward strand.
         * Target without strand (fragments)
         * are also stored in this list.
         */
        std::list<GenomeRegion> target_list_fw ;
        /*!
         * \brief A buffers containing the
         * target mapped on the reverse strand.
         */
        std::list<GenomeRegion> target_list_rv ;

} ;

#endif // CORRELATIONMATRIXCREATOR_HPP
