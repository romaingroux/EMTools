#include<dna_utility.hpp>

#include <string>         // string, to_string()
#include <unordered_map>
#include <stdexcept>      // std::invalid_argument


#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <seqan/seq_io.h>   // seqan::SeqFileIn


int dna::map(char base, bool rev_compl)
{
    static bool init  = false ;
    static std::unordered_map<char,int> hash_map ;
    static std::unordered_map<char,int> hash_map_rev ;
    if(not init)
    {   hash_map['A'] = 0 ;
        hash_map['a'] = 0 ;
        hash_map['C'] = 1 ;
        hash_map['c'] = 1 ;
        hash_map['G'] = 2 ;
        hash_map['g'] = 2 ;
        hash_map['T'] = 3 ;
        hash_map['t'] = 3 ;
        hash_map['N'] = 4 ;
        hash_map['n'] = 4 ;

        hash_map_rev['A'] = hash_map['T'] ;
        hash_map_rev['a'] = hash_map['t'] ;
        hash_map_rev['C'] = hash_map['G'] ;
        hash_map_rev['c'] = hash_map['g'] ;
        hash_map_rev['G'] = hash_map['C'] ;
        hash_map_rev['g'] = hash_map['c'] ;
        hash_map_rev['T'] = hash_map['A'] ;
        hash_map_rev['t'] = hash_map['a'] ;
        hash_map_rev['N'] = hash_map['N'] ;
        hash_map_rev['n'] = hash_map['n'] ;

        init = true ;
    }

    try
    {   if(rev_compl)
        {   return hash_map_rev.at(base) ; }
        else
        {   return hash_map.at(base) ; }
    }
    // key could not be found
    catch(std::out_of_range& e)
    {   char msg[256] ;
        sprintf(msg, "Error! Invalid DNA base : %c", base) ;
        throw std::invalid_argument(msg) ;
    }
}

char dna::map(int base, bool rev_compl)
{
    static bool init  = false ;
    static std::unordered_map<int,char> hash_map ;
    static std::unordered_map<int,char> hash_map_rev ;
    if(not init)
    {   hash_map[0] = 'A' ;
        hash_map[1] = 'C' ;
        hash_map[2] = 'G' ;
        hash_map[3] = 'T' ;
        hash_map[4] = 'N' ;

        hash_map_rev[4] = hash_map[4] ;
        hash_map_rev[3] = hash_map[0] ;
        hash_map_rev[2] = hash_map[1] ;
        hash_map_rev[1] = hash_map[2] ;
        hash_map_rev[0] = hash_map[3] ;

        init = true ;
    }

    try
    {   if(rev_compl)
        {   return hash_map_rev.at(base) ; }
        else
        {   return hash_map.at(base) ; }
    }
    // key could not be found
    catch(std::out_of_range& e)
    {   char msg[256] ;
        sprintf(msg, "Error! Invalid DNA code : %i", base) ;
        throw std::invalid_argument(msg) ;
    }
}

int dna::char_to_int(char c, bool rev_compl)
{   return dna::map(c, rev_compl) ; }

Matrix2D<int> dna::char_to_int(const Matrix2D<int>& matrix)
{
    size_t n_row = matrix.get_nrow() ;
    size_t n_col = matrix.get_ncol() ;

    Matrix2D<int> data_int(n_row, n_col) ;
    for(size_t i=0; i<n_row; i++)
    {   for(size_t j=0; j<n_col; j++)
        {   data_int(i,j) = dna::char_to_int(matrix(i,j)) ; }
    }
    return data_int ;
}

char dna::int_to_char(int n, bool rev_compl)
{   return dna::map(n, rev_compl) ; }

std::vector<double> dna::base_composition(const Matrix2D<int>& sequences, bool both_strands)
{
    double total = 0. ;
    std::vector<double> base_comp(4,0.) ;

    int base_N = dna::map('N') ;

    for(size_t i=0; i<sequences.get_nrow(); i++)
    {   for(size_t j=0; j<sequences.get_ncol(); j++)
        {   // forward strand
            int base = sequences(i,j) ;
            // do not account for N's
            if(base == base_N)
            {   continue ; }
            else
            {   base_comp[base] += 1; 
                total += 1. ;
            }
            // reverse complement strand
            if(both_strands)
            {   // size_t c_hash_rev = dna::hash(c, true) ;
                base_comp[4-base-1] += 1. ;
                total += 1. ;
            }
        }
    }

    // normalize
    for(auto& i : base_comp)
    {   i /= total ; }

    return base_comp ;
}

std::vector<double> dna::base_composition(const Matrix3D<double>& consensus_sequences, bool both_strands)
{
    if(consensus_sequences.get_dim()[2] != 4)
    {   char msg[4096] ;
        sprintf(msg, "Error! consensus sequences 3rd dimension not equal to 4 (%zu)",
                consensus_sequences.get_dim()[2]) ;
        throw std::invalid_argument(msg) ;
    }
    double total = 0. ;
    std::vector<double> base_comp(4,0.) ;

    for(size_t i=0; i<consensus_sequences.get_dim()[0]; i++)
    {   for(size_t j=0; j<consensus_sequences.get_dim()[1]; j++)
        {   for(size_t k=0; k<4; k++)
            {   // forward strand
                {   base_comp[k] += consensus_sequences(i,j,k) ;
                    total += consensus_sequences(i,j,k) ;
                }
                // revers strand
                if(both_strands)
                {   size_t k_comp = 4 - k - 1 ;
                    base_comp[k_comp] += consensus_sequences(i,j,k) ;
                    total += consensus_sequences(i,j,k) ;
                }
            }
        }
    }

    // normalize
    for(auto& i : base_comp)
    {   i /= total ; }

    return base_comp ;
}

std::string dna::consensus_to_kmer(const Matrix3D<double>& consensus_sequences,
                                        size_t row,
                                        size_t from,
                                        size_t to)
{   // dna letter codes
    static int n_code = dna::char_to_int('N') ;

    // some dimensions
    size_t length = to - from ;
    size_t n_dim3 = consensus_sequences.get_dim()[2] ;

    // kmer
    std::string kmer(length, n_code) ;


    for(size_t i_mat=from, i_kmer=0; i_mat<to; i_mat++, i_kmer++)
    {   // get majority base
        int max_j =  5 ;
        double max_n = -1. ;
        for(size_t j=0; j<n_dim3; j++)
        {   if(consensus_sequences(row,i_mat,j) > max_n)
            {   max_n = consensus_sequences(row,i_mat,j) ;
                max_j = j ;
            }
        }
        kmer[i_kmer] = std::to_string(max_j)[0] ;
    }

    return kmer ;
}
