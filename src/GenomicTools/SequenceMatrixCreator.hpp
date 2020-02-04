#ifndef SEQUENCEMATRIXCREATOR_HPP
#define SEQUENCEMATRIXCREATOR_HPP

#include <MatrixCreator.hpp>

#include <seqan/seq_io.h>  // seqan::SeqFileIn
#include <string>
#include <Matrix2D.hpp>

/*!
 * \brief The SequenceMatrixCreator class implements the creation
 * of a DNA sequence matrix with sequences on the rows and individual
 * bases on the columns.
 * The sequences of interest are determined from a BED file and a fasta
 * file containing a genome sequences.
 * The middle of the regions contained in the BED file specify the middle
 * of the sequences which are then extended on both sides until covering
 * a specified <from>/<to> range.
 * The final matrix contains as many sequences as there were regions in
 * the BED file. The first matrix row corresponds to the sequence extracted
 * from the first BED region, and so one.
 * The chromsome/sequence names between the BED and the fasta file must
 * be identical.
 */
class SequenceMatrixCreator : public MatrixCreator
{
    public:
    /*!
     * \brief Constructs an object to build a
     * sequence matrix from a partition.
     * \param bed_file_path the path to the file containing
     * the references.
     * \param fasta_file_path the path to the file containing
     * the sequences.
     * \param from the upstream most relative position
     * to consider around the references. It may
     * be changed to make sure that the central bin
     * is centered on +/- 0.
     * \param to the dowmstream most relative position
     * to consider around the references. It may
     * be changed to make sure that the central bin
     * is centered on +/- 0.
     */
        SequenceMatrixCreator(const std::string& bed_file_path,
                              const std::string& fasta_file_path,
                              int from,
                              int to) ;

        /*!
         * \brief Destructor
         */
        virtual ~SequenceMatrixCreator() ;

        /*!
        * \brief Computes the matrix and returns it.
        * \return the sequence matrix.
        * \throw std::runtime_error if two sequences
        * have the same header in the fasta file or
        * if a sequence/chromosome name present
        * in the bed cannot be found as sequence
        * header in the fasta file.
        */
        virtual Matrix2D<int> create_matrix() override ;

    protected:
        /*!
         * \brief Opens the fasta file.
         * \throw std::runtime_error if the file cannot
         * be open.
         */
        void open_fasta_file() ;

        /*!
         * \brief Closes the fasta file.
         * \throw std::runtime_error if the file cannot
         * be open.
         */
        void close_fasta_file() ;

        /*!
         * \brief Fasta file path.
         */
        std::string fasta_path ;
        /*!
         * \brief An input stream to the
         * fasta file.
         * Use open_fasta_file() to open the stream
         * and close_fasta_file() to close it.
         */
        seqan::SeqFileIn fasta_file ;

} ;


#endif // SEQUENCEMATRIXCREATOR_HPP
