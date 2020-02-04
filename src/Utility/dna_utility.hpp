#ifndef DNA_UTILITY_HPP
#define DNA_UTILITY_HPP

#include <string>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>

namespace dna
{
    /*!
     * \brief Contains the mapping to convert
     * DNA characters to integer codes.
     * Lower and capital characters are accepted.
     * \param base the character of interest.
     * \param rev_compl whether the reverse
     * complement of the character is of interest.
     * \return the corresponding DNA code.
     */
    int map(char base, bool rev_compl=false) ;

    /*!
     * \brief Contains the mapping to convert
     * DNA code to characters.
     * Only capital characters are returned.
     * \param base the code of interest.
     * \param rev_compl whether the reverse
     * complement of the code is of interest.
     * \return the corresponding DNA character.
     */
    char map(int base, bool rev_compl=false) ;

    /*!
     * \brief Converts a DNA character (A, C,
     * G, T) to an integer.
     * \param c the DNA character of interest.
     * \return the character integer code.
     * \throw std::invalid_argument if the
     * given character is not a valid DNA
     * character.
     */
    int char_to_int(char c, bool rev_compl= false) ;

    /*!
     * \brief Converts a DNA character matrix (A, C,
     * G, T) to an integer matrix.
     * The DNA characters are converted using
     * SequenceLayer::char_to_int(char).
     * param file_address the address of the file to load.
     * \return the corresponding int matrix.
     */
    Matrix2D<int> char_to_int(const Matrix2D<int>& matrix) ;

    /*!
     * \brief Converts an int DNA code to
     * a regular DNA character (A, C, G, T).
     * This method is the reverse method of
     * char_to_int(char).
     * \param n the DNA code of interest.
     * \return the DNA character.
     * \throw std::invalid_argument if the
     * given int is not a valid DNA
     * code.
     */
    char int_to_char(int n, bool rev_compl=false) ;

    /*!
     * \brief Computes the base composition of a set of 
     * sequences, in integer format, contained in a matrix.
     * \param sequences a matrix containing the sequences 
     * of interest.
     * \param both_strands also accounts for the reverse 
     * complement of the sequences.
     * \throw std::invalid_argument if a non-supported 
     * character is found in the matrix.
     * \return a vector of 4 values corresponding to the
     * frequencies of A,C,G and T
     * respectively.
     */
    std::vector<double> base_composition(const Matrix2D<int>& sequences,
                                         bool both_strands) ;

    
    /*!
     * \brief Computes the base composition of a set of 
     * consensus sequences, contained in a probability 
     * matrix with the following dimensions:
     * 1st the number of sequences
     * 2nd the length of the sequences
     * 3rd 4 for A,C,G,T
     * The sum over the 1st and 2nd dimension should be 1.
     * The overall sum of the matrix values should be equal
     * to the 1st dimension.
     * \param consensus_sequences a matrix containing the 
     * consensus sequences of interest.
     * \param both_strands also accounts for the reverse 
     * complement of the sequences.
     * \throw std::invalid_argument if the 3rd dimension
     * of the consensus sequence matrix is not equal to 4.
     * \return a vector of 4 values corresponding to the
     * frequencies of A,C,G and T
     * respectively.
     */
    std::vector<double> base_composition(const Matrix3D<double>& consensus_sequences,
                                         bool both_strands) ;

    /*!
     * \brief Constructs a kmer from a piece of consensus
     * sequence piece contained in a probability matrix with
     * the following dimensions:
     * 1st the number of sequences
     * 2nd the length of the sequences
     * 3rd 4 for A,C,G,T
     * The kmer of a consensus sequence is defined as the
     * kmer composed of the most important base at each
     * position.
     * \param consensus_sequences the matrix containing the
     * consensus sequence "kmer".
     * \param row the index of the row containing the
     * consensus sequence piece.
     * \param from the starting position of the consensus
     * sequence piece.
     * \param to the past last ending position of the
     * consensus sequence piece.
     * \return the kmer.
     */
    std::string consensus_to_kmer(const Matrix3D<double>& consensus_sequences,
                                  size_t row,
                                  size_t from,
                                  size_t to) ;
}

#endif // DNA_UTILITY_HPP
