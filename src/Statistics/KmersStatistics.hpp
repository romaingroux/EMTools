#ifndef KMER_STATISTICS_HPP
#define KMER_STATISTICS_HPP

#include <string>
#include <vector>
#include <unordered_map>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>


namespace kmers {

    /*!
     * \brief Creates a table of all kmer found in
     * a set of consensus sequences and return an
     * associated pvalue for each one of them.
     * The pvalue is computed as 1 - p(count, exp)
     * where p is the CDF of a Poisson distribution
     * with <count> being the number of time a given
     * kmer is observed in the data and <exp> its
     * expected number of occurences. The expected
     * number of occurences is computed from the
     * observed di-nucleotide occurences in the given
     * data.
     * \param sequences the sequences in integer
     * format (A:0, C:1, G:2, T:3). It has the
     * following dimensions :
     * 1st the number of sequences
     * 2n the length of the sequences
     * \param k the kmer length.
     * \return a pair of vectors containing the kmers
     * in integer format (A:0, C:1, G:2, T:3) and
     * their pvalues (in the same order).
     */
    std::pair<std::vector<std::string>, std::vector<double>>
    compute_kmer_pvalue(const Matrix2D<int>& sequences,
                        size_t k) ;

    /*!
     * \brief Creates a table of all kmer found in
     * a set of consensus sequences and return an
     * associated pvalue for each one of them.
     * To transform a consensus sequence into a kmer,
     * the major base at each position is considered
     * only.
     * The pvalue is computed as 1 - p(count, exp)
     * where p is the CDF of a Poisson distribution
     * with <count> being the number of time a given
     * kmer is observed in the data and <exp> its
     * expected number of occurences. The expected
     * number of occurences is computed from the
     * observed di-nucleotide occurences in the given
     * data.
     * \param consensus_sequence a matrix containing the
     * consensus sequences as a probability matrix. It has
     * the following dimensions :
     * 1st the number of consensus sequences
     * 2nd the length of the consensus sequences
     * 3rd 4 for A, C, G, T
     * The summing over the 1st and 2nd dimensions should
     * give 1s.
     * \param k the kmer length.
     * \return a pair of vectors containing the kmers
     * in integer format (A:0, C:1, G:2, T:3) and
     * their pvalues (in the same order).
     */
    std::pair<std::vector<std::string>, std::vector<double>>
    compute_kmer_pvalue(const Matrix3D<double>& consensus_sequences,
                        size_t k) ;

    /*!
     * \brief Computes the number of occurences of each kmer
     * appearing in a set of consensus sequences.
     * To transform a consensus sequence into a kmer,
     * the major base at each position is considered
     * only.
     * \param sequences the sequences in integer
     * format (A:0, C:1, G:2, T:3). It has the
     * following dimensions :
     * 1st the number of sequences
     * 2n the length of the sequences
     * \param k the kmer length.
     * \return a map with the kmer found as keys and their
     * number of occurences as value.
     */
    std::unordered_map<std::string, int>
    compute_kmer_count(const Matrix2D<int>& sequences,
                       size_t k) ;

    /*!
     * \brief Computes the number of occurences of each kmer
     * appearing in a set of consensus sequences.
     * To transform a consensus sequence into a kmer,
     * the major base at each position is considered
     * only.
     * \param consensus_sequence a matrix containing the
     * consensus sequences as a probability matrix. It has
     * the following dimensions :
     * 1st the number of consensus sequences
     * 2nd the length of the consensus sequences
     * 3rd 4 for A, C, G, T
     * The summing over the 1st and 2nd dimensions should
     * give 1s.
     * \param k the kmer length.
     * \return a map with the kmer found as keys and their
     * number of occurences as value.
     */
    std::unordered_map<std::string, int>
    compute_kmer_count(const Matrix3D<double>& consensus_sequence,
                       size_t k) ;

    /*!
     * \brief Compute the pvalue of a set of kmers given their
     * number of occurences and their expected number of
     * occurences.
     * \param counts the number of occurences of each kmer.
     * \param expected the expected number of occurences of
     * each kmer.
     * \return the pvalue of each kmer.
     */
    std::vector<double> compute_pvalues(const std::vector<int>& counts,
                                        const std::vector<double>& expected) ;

    /*!
     * \brief Computes the number of expected occurences of a set
     * of kmers given the mononucleotide and dinucleotide composition
     * of the kmers.
     * \param k the length of the kmers.
     * \param kmers the kmers in integer format (A:0, C:1, G:2, T:3)
     * \param comp1 a matrix with the number of counts of each
     * base at each position in the kmers. It has the following
     * dimensions :
     * 1st dim k
     * 2nd dim 4 for A, C, G, T
     * \param comp2 a matrix containing the number of counts of
     * each pair of dinucleotide in the kmers. It has the
     * following dimensions :
     * 1st 4 for A, C, G, T (1st base of the dinucleotide)
     * 2nd k-1
     * 3rd 4 for A, C, G, T (2n base of the dinucleotide)
     * \return the expected number of occurences for each
     * kmer.
     */
    std::vector<double> compute_exp_values(size_t k,
                                          const std::vector<std::string>& kmers,
                                          const Matrix2D<int>& comp1,
                                          const Matrix3D<int>& comp2) ;

    /*!
     * \brief Computes the number of occurences of each dinucleotide
     * pair at each position of a given set of kmers.
     * \param k the length of the kmers.
     * \param kmers the kmers in integer format (A:0, C:1, G:2, T:3).
     * \param counts the number of occurences of each kmer.
     * \return a matrix containing the number of counts of
     * each pair of dinucleotide in the kmers. It has the
     * following dimensions :
     * 1st 4 for A, C, G, T (1st base of the dinucleotide)
     * 2nd k-1
     * 3rd 4 for A, C, G, T (2n base of the dinucleotide)
     */
    Matrix3D<int> compute_dinucleotide(size_t k,
                                       const std::vector<std::string>& kmers,
                                       const std::vector<int>& counts) ;

    /*!
     * \brief Computes the number of occurences of each nucleotide
     * at each position of a given set of kmers.
     * \param k the length of the kmers.
     * \param kmers the kmers in integer format (A:0, C:1, G:2, T:3).
     * \param counts the number of occurences of each kmer
     * \return a matrix with the number of counts of each
     * base at each position in the kmers. It has the following
     * dimensions :
     * 1st dim k
     * 2nd dim 4 for A, C, G, T
     */
    Matrix2D<int> compute_mononucleotide(size_t k,
                            const std::vector<std::string>& kmers,
                            const std::vector<int>& counts) ;

}

#endif // KMER_STATISTICS_HPP
