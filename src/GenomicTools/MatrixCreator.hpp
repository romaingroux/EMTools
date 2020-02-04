#ifndef MATRIXCREATOR_HPP
#define MATRIXCREATOR_HPP

#include <string>

#include <seqan/bed_io.h>  // BedFileIn, BedRecord

#include <GenomeRegion.hpp>
#include <Matrix2D.hpp>


/*!
 * \brief The MatrixCreator class is a base class
 * to be derived by classes that are dedicated to
 * construct data matrices which rows contains
 * a signal at different positions (columns) in
 * this given region.
 */
class MatrixCreator
{   public:
        /*!
         * \brief Returns the central position of a bed region.
         * \param bed_line the region of interest.
         * \return the position of the center.
         */
        static int get_center_pos(const seqan::BedRecord<seqan::Bed3>& bed_line) ;

    public:
        /*!
         * \brief Constructs an object.
         * \param bed_file_path the path to the bed file
         * containing the coordinates of the regions of
         * interest.
         * \param from the downstream most position
         * to consider, relative to a set of genomic
         * positions.
         * \param to the upstream most position to
         * consider, relative to a set of genomic
         * positions
         */
        MatrixCreator(const std::string& bed_file_path,
                      int from,
                      int to) ;
        /*!
         * Destructor.
         */
        virtual ~MatrixCreator() ;

        /*!
         * \brief Creates and return the count matrix.
         * \return the count matrix.
         */
        virtual Matrix2D<int> create_matrix() = 0 ;

    protected:

        /*!
         * \brief Opens the bed file.
         * \throw std::runtime_error if the file cannot
         * be open.
         */
        void open_bed_file() ;

        /*!
         * \brief Closes the bed file.
         * Does nothing if already closed.
         */
        void close_bed_file() ;

        /*!
         * \brief Bed file path.
         */
        std::string bed_path ;
        /*!
         * \brief An input stream to the
         * bed file.
         * Use open_bed_file() to open the stream
         * and close_bed_file() to close it.
         */
        seqan::BedFileIn bed_file ;

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
         * \brief A matrix containing the data
         * found at each position around each reference.
         * This is the data structure to fill.
         */
        Matrix2D<int> matrix ;
} ;


#endif // MATRIXCREATOR_HPP


