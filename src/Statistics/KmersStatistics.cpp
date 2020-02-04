#include <KmersStatistics.hpp>

#include <string>
#include <vector>
#include <unordered_map>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Statistics.hpp>    // poisson_cdf()
#include <dna_utility.hpp>   // dna::consensus_to_kmer()

Matrix2D<int> kmers::compute_mononucleotide(size_t k,
                                            const std::vector<std::string>& kmers,
                                            const std::vector<int>& counts)
{
    Matrix2D<int> comp(k, 4, 0) ;

    for(size_t i=0; i<kmers.size(); i++)
    {
        std::string kmer = kmers[i] ;
        int count        = counts[i] ;

        for (size_t j=0; j<k ; j++)
        {   int base = int(kmer[j]) - int('0') ;
            comp(j, base) += count ;
        }
    }
    return(comp);
}

Matrix3D<int> kmers::compute_dinucleotide(size_t k,
                                          const std::vector<std::string>& kmers,
                                          const std::vector<int>& counts)
{
    Matrix3D<int> comp2(4, k-1, 4) ;

    for(size_t i=0; i<kmers.size(); i++)
    {   std::string kmer = kmers[i] ;
        int count        = counts[i] ;
        for (size_t j=0; j<k-1; j++)
        {   int base1 = int(kmer[j])   - int('0') ;
            int base2 = int(kmer[j+1]) - int('0') ;
            comp2(base1, j, base2) += count ;
        }
    }
    return(comp2);
}

std::vector<double> kmers::compute_exp_values(size_t k,
                                              const std::vector<std::string>& kmers,
                                              const Matrix2D<int>& comp1,
                                              const Matrix3D<int>& comp2)
{   std::vector<double> expected(kmers.size()) ;

    for(size_t i=0; i<kmers.size(); i++)
    {   std::string kmer = kmers[i] ;
        double exp       = comp1(0, int(kmer[0]) - int('0')) ;

        for(size_t j=0; j<k-1; j++)
        {   int base1 = int(kmer[j])   - int('0') ;
            int base2 = int(kmer[j+1]) - int('0') ;
            exp      *= ((double)comp2(base1, j, base2) /
                         (double)comp1(j+1, base2)) ;
        }

        expected[i] = exp ;
    }
    return expected;
}


std::vector<double> kmers::compute_pvalues(const std::vector<int>& counts,
                                           const std::vector<double>& expected)
{
    std::vector<double> pvalues(counts.size()) ;

    for(size_t i=0; i<counts.size(); i++)
    {
        pvalues[i] = poisson_cdf((double)counts[i],
                                 expected[i],
                                 false) ;
    }
    return pvalues ;
}

std::unordered_map<std::string, int>
kmers::compute_kmer_count(const Matrix2D<int>& sequences,
                           size_t length)
{   std::unordered_map<std::string, int> kmer_count ;

    size_t n_row   = sequences.get_dim()[0] ;
    size_t n_shift = sequences.get_dim()[1] - length + 1 ;

    for(size_t i=0; i<n_row; i++)
    {   for(size_t s=0; s<n_shift; s++)
        {   // get kmer
            std::string kmer(length, '0') ;
            bool n_found = false ;
            for(size_t j_seq=s, j_kmer=0; j_seq<s+length; j_seq++, j_kmer++)
            {   // check for N's
                if(sequences(i, j_seq) == dna::char_to_int('N'))
                {   n_found = true ;
                    break ;
                }
                kmer[j_kmer] = std::to_string(sequences(i, j_seq))[0] ;
            }

            // found N -> skip this one
            if(n_found)
            {   continue ; }

            // insert
            auto iter = kmer_count.find(kmer) ;
            if(iter == kmer_count.end())
            {   kmer_count.emplace(kmer, 1) ; }
            // update
            else
            {   iter->second += 1 ; }
        }
    }
    return kmer_count ;
}

std::unordered_map<std::string, int>
kmers::compute_kmer_count(const Matrix3D<double>& consensus_sequence,
                           size_t length)
{   std::unordered_map<std::string, int> kmer_count ;

    size_t n_row   = consensus_sequence.get_dim()[0] ;
    size_t n_shift = consensus_sequence.get_dim()[1] - length + 1 ;

    for(size_t i=0; i<n_row; i++)
    {   for(size_t s=0; s<n_shift; s++)
        {   // get kmer
            std::string kmer = dna::consensus_to_kmer(consensus_sequence,
                                                      i,
                                                      s,
                                                      s+length) ;
            // insert
            auto iter = kmer_count.find(kmer) ;
            if(iter == kmer_count.end())
            {   kmer_count.emplace(kmer, 1) ; }
            // update
            else
            {   iter->second += 1 ; }
        }
    }
    return kmer_count ;
}

std::pair<std::vector<std::string>, std::vector<double>>
kmers::compute_kmer_pvalue(const Matrix2D<int>& sequences,
                           size_t k)
{

    // get kmer count
    std::unordered_map<std::string, int> kmer_count =
            kmers::compute_kmer_count(sequences, k) ;

    // turn to vectors
    std::vector<std::string> kmers(kmer_count.size()) ;
    size_t i=0 ;
    for(const auto& iter : kmer_count)
    {   kmers[i] = iter.first ;
        i++ ;
    }
    std::sort(kmers.begin(), kmers.end()) ;
    std::vector<int> counts(kmer_count.size()) ;
    i=0 ;
    for(const auto& kmer : kmers)
    {   counts[i] = kmer_count[kmer] ;
        i++ ;
    }
    kmer_count.clear() ;

    // get mononucleotide composition in kmer
    // this is a probability matrix modeling the kmer
    Matrix2D<int> comp1 = kmers::compute_mononucleotide(k, kmers, counts) ;

    // get dinucleotide composition in kmer
    // this is a probability matrix modeling the kmer
    Matrix3D<int> comp2 = kmers::compute_dinucleotide(k, kmers, counts) ;

    // get expected number of occurence for each kmer
    std::vector<double> expected = kmers::compute_exp_values(k,
                                                             kmers,
                                                             comp1,
                                                             comp2) ;

    // compute the pvalue associated to each kmer
    std::vector<double> pvalues = kmers::compute_pvalues(counts,
                                                         expected) ;

    return std::make_pair(kmers, pvalues) ;
}

std::pair<std::vector<std::string>, std::vector<double>>
kmers::compute_kmer_pvalue(const Matrix3D<double>& consensus_sequences,
                    size_t k)
{

    // get kmer count
    std::unordered_map<std::string, int> kmer_count =
            kmers::compute_kmer_count(consensus_sequences, k) ;

    /*
    consensus_sequences.get_dim() ;
    std::unordered_map<std::string, int> kmer_count ;
    // ANN
    kmer_count["000"] = 1  ;
    kmer_count["001"] = 2  ;
    kmer_count["002"] = 3  ;
    kmer_count["003"] = 4  ;
    kmer_count["010"] = 5  ;
    kmer_count["011"] = 6  ;
    kmer_count["012"] = 7  ;
    kmer_count["013"] = 8  ;
    kmer_count["020"] = 9  ;
    kmer_count["021"] = 10 ;
    kmer_count["022"] = 11 ;
    kmer_count["023"] = 12 ;
    kmer_count["030"] = 13 ;
    kmer_count["031"] = 14 ;
    kmer_count["032"] = 15 ;
    kmer_count["033"] = 16 ;
    // CNN
    kmer_count["100"] = 1  ;
    kmer_count["101"] = 2  ;
    kmer_count["102"] = 3  ;
    kmer_count["103"] = 4  ;
    kmer_count["110"] = 5  ;
    kmer_count["111"] = 6  ;
    kmer_count["112"] = 7  ;
    kmer_count["113"] = 8  ;
    kmer_count["120"] = 9  ;
    kmer_count["121"] = 10 ;
    kmer_count["122"] = 11 ;
    kmer_count["123"] = 12 ;
    kmer_count["130"] = 13 ;
    kmer_count["131"] = 14 ;
    kmer_count["132"] = 15 ;
    kmer_count["133"] = 16 ;
    // GNN
    kmer_count["200"] = 1  ;
    kmer_count["201"] = 2  ;
    kmer_count["202"] = 3  ;
    kmer_count["203"] = 4  ;
    kmer_count["210"] = 5  ;
    kmer_count["211"] = 6  ;
    kmer_count["212"] = 7  ;
    kmer_count["213"] = 8  ;
    kmer_count["220"] = 9  ;
    kmer_count["221"] = 10 ;
    kmer_count["222"] = 11 ;
    kmer_count["223"] = 12 ;
    kmer_count["230"] = 13 ;
    kmer_count["231"] = 14 ;
    kmer_count["232"] = 15 ;
    kmer_count["233"] = 16 ;
    // TNN
    kmer_count["300"] = 1  ;
    kmer_count["301"] = 2  ;
    kmer_count["302"] = 3  ;
    kmer_count["303"] = 4  ;
    kmer_count["310"] = 5  ;
    kmer_count["311"] = 6  ;
    kmer_count["312"] = 7  ;
    kmer_count["313"] = 8  ;
    kmer_count["320"] = 9  ;
    kmer_count["321"] = 10 ;
    kmer_count["322"] = 11 ;
    kmer_count["323"] = 12 ;
    kmer_count["330"] = 13 ;
    kmer_count["331"] = 14 ;
    kmer_count["332"] = 15 ;
    kmer_count["333"] = 16 ;
    */

    // turn to vectors
    std::vector<std::string> kmers(kmer_count.size()) ;
    size_t i=0 ;
    for(const auto& iter : kmer_count)
    {   kmers[i] = iter.first ;
        i++ ;
    }
    std::sort(kmers.begin(), kmers.end()) ;
    std::vector<int> counts(kmer_count.size()) ;
    i=0 ;
    for(const auto& kmer : kmers)
    {   counts[i] = kmer_count[kmer] ;
        i++ ;
    }
    kmer_count.clear() ;

    // get mononucleotide composition in kmer
    // this is a probability matrix modeling the kmer
    Matrix2D<int> comp1 = kmers::compute_mononucleotide(k, kmers, counts) ;

    // get dinucleotide composition in kmer
    // this is a probability matrix modeling the kmer
    Matrix3D<int> comp2 = kmers::compute_dinucleotide(k, kmers, counts) ;

    // get expected number of occurence for each kmer
    std::vector<double> expected = kmers::compute_exp_values(k,
                                                             kmers,
                                                             comp1,
                                                             comp2) ;

    // compute the pvalue associated to each kmer
    std::vector<double> pvalues = kmers::compute_pvalues(counts,
                                                         expected) ;

    return std::make_pair(kmers, pvalues) ;
}
