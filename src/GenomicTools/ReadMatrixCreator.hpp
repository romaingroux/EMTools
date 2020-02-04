#ifndef READMATRIXCREATOR_HPP
#define READMATRIXCREATOR_HPP

#include <MatrixCreator.hpp>

#include <unordered_map>
#include <utility>         // std::pair, std::make_pair()

#include <seqan/bed_io.h>  // BedFileIn
#include <seqan/bam_io.h>  // BamFileIn, BamAlignmentRecord
#include <GenomeRegion.hpp>

/*!
 * \brief The ReadMatrixCreator class provides an interface for classes
 * that will be dedicated to the creation of read density matrices.
 * It provides an interface to open/close BAM and BAI files and to
 * handle sequencing reads and sequencing fragments.
 */
class ReadMatrixCreator : public MatrixCreator
{
    public:
        /*!
         * \brief A list of values indicating how the data
         * should be handled when counting the number of
         * fragments mapped in a given bin.
         *
         * FRAGMENT : all positions within a fragment are
         * accounted for and attributed to the
         * corresponding bins :
         *       bin1    bin2
         * ----|-------|-------|------------> genome
         *   -------  -------                 fragments
         *   --> <--  --> <--                 pair of reads
         *     |||||   ||||||                 scoring positions
         * bin1 gets a score of 5 and bin2 a
         * score of 6.
         *
         * FRAGMENT_CENTER : only the central position
         * within a fragment is accounted for and
         * attributed to the corresponding bin :
         * *       bin1    bin2
         * ----|-------|-------|------------> genome
         *   -------  -------                 fragments
         *   --> <--  --> <--                 pair of reads
         *      |        |                    scoring positions
         * bin1 gets a score of 1 and bin2 also.
         *
         * READ : all positions within a read are
         * accounted for and attributed to the
         * corresponding bins :
         *       bin1    bin2
         * ----|-------|-------|------------> genome
         *   -------    -------               fragments
         *   --> <--    --> <--               reads
         *     | |||    ||| |||               scoring positions
         * bin1 gets a score of 4 and bin2 a
         * score of 6.
         *
         * READ_ATAC : only the shifted start
         * of the reads are used. Additionally, the
         * start position is shifted by +4bp(towards
         * the right) for reads on the + strand and
         * -5bp for reads on the - strand (towards the
         * left). These positions indicate the insertion
         * position in ATAC-seq data.
         *       bin1    bin2
         * ----|-------|-------|------------> genome
         *   -------    -------               fragments
         *   --> <--    --> <--               reads
         *         |    |     |               scoring positions
         * bin1 gets a score of 1 and bin2 a
         * score of 2.
         */
        enum methods {FRAGMENT=0,
                      FRAGMENT_CENTER,
                      READ,
                      READ_ATAC} ;

    public:

        /*!
         * \brief Computes which bins (from a contiguous
         * range of bins) are overlapped by a given target
         * and returns two indices corresponding to :
         * i) the index of the 1st bin overlapped by the
         * target
         * ii) the index of the past last bin overlapepd
         * by the target.
         * If the target does not overlapp any bin (it is
         * located upstream the 1st bin, downstream the
         * last bin or on a different chromosome), the
         * index pair 0,0  is returned.
         * Thus, in any case, a loop of the type
         * for(i=first,i<second,i++) can always be used
         * on the result.
         * \param target t
         * \param bins
         * \return
         */
        static std::pair<int, int> get_bin_indices(const GenomeRegion& target,
                                                   const std::vector<GenomeRegion>& bins) ;

        /*!
         * \brief Checks that the read is i) is mapped
         * , ii) passes QC and iii) is not a duplicate,
         * based on the flag value.
         * \param read the read of interest.
         * \return whether the read passes the above tests.
         */
        bool is_good_read(const seqan::BamAlignmentRecord& read) ;

        /*!
         * \brief Checks that the read is i) a good read, ii)
         * a paired read, iii) proplery aligned,  iv) the 1st
         * of the pair based on the flag values and that
         * v) they forms a proper fragment with its mate mate
         * (both read should point toward one other).
         * \param read the read of interest.
         * \return whether the read and its mate form a proper
         * fragment.
         */
        bool is_good_pair(const seqan::BamAlignmentRecord& read) ;

    public:

        ReadMatrixCreator() = delete ;

        /*!
         * \brief Constructs an object to create
         * a genomic count matrix.
         * \param bed_file_path the path to the file containing
         * the references.
         * \param bam_file_path the path to the file containing
         * the targets.
         * \param bai_file_path the path to index file of the bam
         * file containing the targets.
         * \param from the downstream most position
         * to consider, relative to a set of genomic
         * positions.
         * \param to the upstream most position to
         * consider, relative to a set of genomic
         * positions
         * \param bin_size the size of the bins in
         * which the regions encompassing the set
         * of genomic positions will be broken
         * into.
         * \param method how the sequenced fragments
         * should be consider when assigning counts
         * to the bins.
         */
        ReadMatrixCreator(const std::string& bed_file_path,
                          const std::string& bam_file_path,
                          const std::string& bai_file_path,
                          int from,
                          int to,
                          int bin_size,
                          ReadMatrixCreator::methods method) ;

        /*!
         * Destructor.
         */
        virtual ~ReadMatrixCreator() ;

    protected:

        /*!
         * \brief Binarize the given range [from,to] into
         * equal sized bins having the specified size.
         * The bin coordinates are stored in bin_coord as
         * pairs of [start,end) coordinates. One bin is
         * centered on +/- 0.
         * This piece of code has been taken from ChIP-seq
         * source code on SourceForge (from ChIPExtract
         * program).
         */
        void compute_relative_bin_coord() ;

        /*!
         * \brief Checks whether a record has a valid chromosome,
         * that is whether this chromosome has been found in the
         * bed file has well.
         * \param record a record from the bam file.
         * \return whether the record chromosome is valid.
         */
        bool is_valid_chromosome(const seqan::BamAlignmentRecord& record) ;

        /*!
         * \brief Opens the bam file.
         * \throw std::runtime_error if the file cannot
         * be open.
         */
        void open_bam_file() ;

        /*!
         * \brief Opens the bam index file.
         * \throw std::runtime_error if the file cannot
         * be open.
         */
        void open_bai_file() ;

        /*!
         * \brief Closes the bam file.
         * Does nothing if already closed.
         */
        void close_bam_file() ;

        /*!
         * \brief The bin size.
         */
        int bin_size ;
        /*!
         * \brief How to consider the sequenced fragments when computing
         * the bin values.
         */
        ReadMatrixCreator::methods method ;
        /*!
         * \brief The relative bin coordinates, compared to a given
         * position. Each bin has a pair [from,to) where <from> is the
         * 1st position within the bin and <to> is the 1st position
         * after the bin. One bin is centered on +/- 0.
         */
        std::vector<std::pair<int,int>> relative_bin_coord ;

        /*!
         * \brief Bam file path.
         */
        std::string bam_path ;
        /*!
         * \brief Bam index file path.
         */
        std::string bai_path ;

        /*!
         * \brief An input stream to the
         * bam file.
         * Use open_bam_file() to open the stream
         * and close_bam_file() to close it.
         */
        seqan::BamFileIn bam_file;
        /*!
         * \brief An input stream to the
         * bam index file.
         * Use open_bai_file() to open the stream
         * and close_bai_file() to close it.
         */
        seqan::BamIndex<seqan::Bai> bai_file ;
        /*!
         * \brief A map containing the valid chromsome
         * names as keys (as find in the bed file) and
         * their indices (as found in the BAM header)
         * as values.
         */
        std::unordered_map<std::string, int> chrom_map_names ;
        /*!
         * \brief A vector containing containing,
         * for each reference, the coordinates of
         * the genomic region covered by the bins.
         */
        std::vector<std::vector<GenomeRegion>> matrix_bins ;
} ;




#endif // READMATRIXCREATOR_HPP
