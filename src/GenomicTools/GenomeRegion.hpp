#ifndef GENOMEREGION_HPP
#define GENOMEREGION_HPP

#include <seqan/bam_io.h>
#include <iostream>
#include <string>

/*!
 * \brief The GenomeRegion class models a segment along
 * a genome. It encompasses one or more base pairs.
 * It is characterized by a pair of start/end positions
 * defining an open range [start,end) with start being the
 * 1st position within the region and end the 1st position
 * outside. All coordinates are expected to be 0-based.
 * A genome region is unoriented - no strand information - and
 * as such, the starting coordinate should always be smaller
 * than the end coordinate.
 */
class GenomeRegion
{
    public:
        /*!
         * \brief Returns an object corresponding to
         * the sequenced read.
         * The read is assumed to be of proper quality. No check
         * is performed.
         * \param record an alignment read from a BAM file.
         * \param bam_file an open stream to the corresponding
         * BAM file.
         * \return a region covering the sequenced read alignment.
         * \throw std::invalid_argument if a position
         * is negative of if end is smaller or equal
         * to start.
         */
        static GenomeRegion constructRead(const seqan::BamAlignmentRecord& record,
                                          const seqan::BamFileIn& bam_file) ;

        /*!
         * \brief Returns an object corresponding to the
         * transposition site for ATAC-seq data. It corresponds
         * to the primary read start position shifted by +4bp (to the
         * right) for reads on the + strand and -5bp (to the left)
         * for reads on the - strand.
         * The read is assumed to be of proper quality. No check
         * is performed.
         * \\param record an alignment read from a BAM file.
         * \param bam_file an open stream to the corresponding
         * BAM file.
         * \return a region corresponding to the transposition
         * site for ATAC-seq reads.
         */
        static GenomeRegion constructReadATAC(const seqan::BamAlignmentRecord& record,
                                              const seqan::BamFileIn& bam_file) ;

        /*!
         * \brief Returns an object corresponding to the
         * starting edge of the read.
         * The read is assumed to be of proper quality. No check
         * is performed.
         * \\param record an alignment read from a BAM file.
         * \param bam_file an open stream to the corresponding
         * BAM file.
         * \return a region corresponding to the edgeof the
         * read.
         */
        static GenomeRegion constructReadEdge(const seqan::BamAlignmentRecord& record,
                                              const seqan::BamFileIn& bam_file) ;

        /*!
         * \brief Returns an object corresponding to
         * the fragment contained by the two sequenced
         * reads of the pair (the BAM file should contain
         * paired-end reads).
         * The read is assumed to be properly paired with
         * another reads and no check is performed.
         * \param record an alignment read from a BAM file.
         * \param bam_file an open stream to the corresponding
         * BAM file.
         * \return a region covering the sequenced fragment
         * alignment.
         * \throw std::invalid_argument if a position
         * is negative of if end is smaller or equal
         * to start.
         */
        static GenomeRegion constructFragment(const seqan::BamAlignmentRecord& record,
                                              const seqan::BamFileIn& bam_file) ;

        /*!
         * \brief Returns an object corresponding to
         * the central position of the fragment contained
         * by the two sequenced reads of the pair (the BAM
         * file should contain paired-end reads).
         * The read is assumed to be properly paired with
         * another reads and no check is performed.
         * \param record an alignment read from a BAM file.
         * \param bam_file an open stream to the corresponding
         * BAM file.
         * \return a region covering the central position of the
         * sequenced fragment alignment.
         * \throw std::invalid_argument if a position
         * is negative of if end is smaller or equal
         * to start.
         */
        static GenomeRegion constructFragmentCenter(const seqan::BamAlignmentRecord& record,
                                                    const seqan::BamFileIn& bam_file) ;
    public:
        /*!
         * Constructs an empty object.
         */
        GenomeRegion() = default ;

        /*!
         * \brief Copy constructor.
         * \param other the other genomic region
         */
        GenomeRegion(const GenomeRegion& other) ;

        /*!
         * \brief Creates a genomic region at the
         * given position.
         * \param chromosome the name of the
         * chromosome on which the region is.
         * \param start the 1st position, 0 based,
         * in the region.
         * \param end the 1st position, 0 based,
         * after the region.
         * \param chromosome_idx the index of the
         * chromosome on which the region is.
         * \throw std::invalid_argument if a position
         * is negative of if end is smaller or equal
         * to start.
         */
        GenomeRegion(const std::string& chromosome,
                     int chromosome_idx,
                     int start,
                     int end) ;

        /*!
         * \brief Returns whether both regions
         * overlap by at least one position.
         * \param other the second region.
         * \return whether there is an overlap.
         */
        bool overlap(const GenomeRegion& other) const ;

        /*!
         * \brief Returns the length of the overlap
         * between both regions.
         * \param other the second region.
         * \return the length of the overlap
         * between both regions.
         */
        int overlap_len(const GenomeRegion& other) const ;

        /*!
         * \brief Assignment operator.
         * \param rhs the other region.
         * \return a reference to the current object.
         */
        GenomeRegion& operator = (const GenomeRegion& rhs) ;

        /*!
         * \brief Move assignment operator.
         * \param rhs the other region.
         * \return a reference to the current object.
         */
        GenomeRegion& operator = (GenomeRegion&& rhs) ;

        /*!
         * \brief Checks equality.
         * \param rhs the other object.
         * \return whether both objects are the
         * same.
         */
        bool operator == (const GenomeRegion& rhs) const ;

        /*!
         * \brief Overlap operator.
         * Returns whether both regions
         * overlap by at least one position.
         * \param other the other region.
         * \return whether there is an overlap.
         */
        bool operator | (const GenomeRegion& rhs) const ;

        /*!
         * \brief Is upstream operator.
         * Checks if the current object is
         * located stricly upstream (no overlap)
         * of the other region.
         * The chromosome index values are also
         * tested. A smaller chromosome value is
         * interpreted as upstream.
         * \param rhs the other region.
         * \return if the current region is upstream,
         * without overlapping the other region.
         */
        bool operator < (const GenomeRegion& rhs) const ;

        /*!
         * \brief Is downstream operator.
         * Checks if the current object is
         * located stricly downstream (no overlap)
         * of the other region.
         * The chromosome index values are also
         * tested. A smaller chromosome value is
         * interpreted as upstream.
         * \param rhs the other region.
         * \return if the current region is downstream,
         * without overlapping the other region.
         */
        bool operator > (const GenomeRegion& rhs) const ;

        /*!
         * \brief The chromosome name on which the
         * region is. This field is present for
         * representations aesthetic only and is
         * never used otherwise.
         */
        std::string chromosome ;
        /*!
         * \brief The index of the chromosome on which the
         * region is. Indexes are used for sorting per
         * position. Chromosomes with smaller indexes
         * are considered upstream to chromosomes with
         * higher indexes.
         */
        int chromosome_idx ;
        /*!
         * \brief The 1st position, 0 based, within the region.
         */
        int start ;
        /*!
         * \brief The 1st position, 0 based, after the region.
         */
        int end ;

        /*!
         * \brief The number of positions within the region.
         */
        int length ;
} ;

/*!
 * \brief Sends a representation of the region to the given
 * stream.
 * \param stream the stream of interest.
 * \param region the region of interest.
 * \return a reference to the stream.
 */
std::ostream& operator << (std::ostream& stream, const GenomeRegion& region) ;

#endif // GENOMEREGION_HPP
