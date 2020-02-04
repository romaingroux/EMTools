EMTools
========

This repository contains C++ classe and wrapper applications that allow to partition genomic regions using different probabilistic algorithms and to manipulate the results. The partitioning algorithms allow to group regions based on i) the densities of reads (for a given experiment) that mapped on the regions of interest, ii) the DNA sequences of the regions of interest, iii) one or more read densities and (optionally) the DNA sequences of the regions of interest.

# 1 Repository content overview

This repository is composed of the following folders and files :

* images/ contains the images displayed below
* res/ contains some necessary material the R code displaying sequence models
* src/ contains the C++ source code
	* Applications/ contains the implementations of the autonomous applications described below.
	* Clustering/ contains classes implementating the different classifications procedures described below and classes computing the class models.
	* GenomicTools/ contains classes implementing different procedures to create read density or DNA sequence matrices. 
	* GUI/ contains classes to display informations. 
	* Matrix/ contains matrix template classes.
	* Parallel/ contains a class implementing a thread pool.
	* Random contains functions and classes to generate random numbers.
	* Statistics contains different statistical functions.
	* Unittests/ contains unit tests
	* Utility/ contains miscellaneous auxilary functions for DNA, matrix and sorting operations.
* scripts/ contains shell scripts to install the required libraries
* build.sh a shell script to compile everyhing
* LICENSE.md the license file
* README.md the document containing this text

# 2 Data

Each of the partitioning algorithm takes as input a matrix that contains, on each row, a signal over a region. Two main types of data matrices exist : i) read density matrices and ii) sequence matrices

![](images/matrix.png?raw=true)
**Data matrix creation process**
Read density matrices are created using a BED, a BAM and a BAI file. The central position of each region present in the BED is considered (the dashed lines) and are aligned in the final matrix. Each BED region is represented by a row in the matrix (the 1st region is on the top row, the last region is on the bottom row). For each region, a bin - which size, in bp, is specified by the user - is centered on this position. Non-overlapping bins are then constructed on each side (represented by the color rectangles) and the number of reads/fragments/fragment centers - as defined by the user - mapping in each bin are counted and the values are written in the corresponding matrix cells.
DNA sequence matrices are created using a BED and a FASTA file. As for read density matrices, the central positions are aligned in the matrix, with one sequence per row. The sequences in these regions extracted and each nucleotide is written in a single cell. Thus the bin size is forced to 1bp. The characters are encoded as integers as follows : A->0, C->1, G->2, T->3, N->4.

The regions from which the signal - read density or sequence - will be measured are determined using a BED file and a from/to range. First, the middle position of all regions listed in the BED file will be computed. These positions will be the central positions of each region in the matrix. The central bin of each region is constructed around these center positions. The center position is always in the center of the central bins. The surrounding bins are then constructed on each side of the central bins until covering the entire range covered by the from/to parameters. Because the central bin is centered on the central positions determined from the BED file, the boundaries of each regions in the matrix may not exactly fit the from/to range.

A special care should be taken to ensure that i) the chromosome names are strictly identical between the BED, the BAM and the FASTA files and that ii) the BED and BAM are sorted by chromosomes and positions.

## 2.1 Read densities

Creating a partition of regions based on their sequencing profiles (generated using any sequencing assay such as ChIP-seq, DNase-seq, MNase-seq, ATAC-seq, etc), also called read densities, necessitates to construct a 2D matrix of dimensions NxL. This matrix contains N regions (over each row) sub-divided into L bins. Each bin covers a constant number of base-pair, for instance 10bp. The bins strictly follow each other and do not overlap each other.

Read density matrices can be generated using the CorrelationMatrixCreator class in GenomicTools. The following snippet creates a density matrix with regions spanning 2kb, divided in bins of 10bp. The "read" parameter indicates that only the 1st position of the reads is the signal of interest. The other parameters are "fragment", "fragment_center" and "atac_read" that consider the entire fragments (for paired-end sequencing data), the fragment central position (for paired-end sequencing data) and 1st position of reads offseted (+4bp for reads mapping on the +strand and -5bp for reads mapping on the -strand) as the signal of interest respectively.

```
#include <iostream>
#include <string>
#include "Matrix/Matrix2D.hpp"
#include "GenomicTools/CorrelationMatrixCreator.hpp"

CorrelationMatrixCreator mc("/some/bed/file.bed",
                            "/some/bam/file.bam",
                            "/some/bai/file.bai",
                            -1000,
                             1000,
                               10,
                            "read") ;
Matrix2D<int> m = mc.create_matrix() ;
std::cout << m << std::endl ;
```
The same can be achieved by using the standalone application CorrelationMatrixCreator which is a wrapper around an instance of CorrelationMatrixCreator.

```
CorrelationMatrixCreator --bed /some/bed/file.bed --bam /some/bam/file.bam --bai /some/bai/file.bai --from -1000 --to 1000 --binSize 10 --method read
```

For help, enter :

```
CorrelationMatrixCreator -h
```


## 2.2 DNA sequences

To create a partition of regions based on their DNA sequences necessitates to construct a 2D matrix of dimensions NxL. This matrix contains N region sequences (over each row) of length L. Contrarily to the read density matrix, the bin size is 1bp. That is, each cell will contain a single nucleotide.

DNA characters are encoded in integers as follows : A->0, C->1, G->2, T->3, N->4.

DNA sequence matrices can be generated using the SequenceMatrixCreator class in GenomicTools. The following snippet creates a sequence matrix with sequences of 200bp.

```
#include <iostream>
#include <string>
#include "Matrix/Matrix2D.hpp"
#include "GenomicTools/SequenceMatrixCreator.hpp"
SequenceMatrixCreator mc("/some/bed/file.bed",
                         "/some/fasta/file.fasta",
                          -100,
                           100) ;
Matrix2D<int> m = mc.create_matrix() ;
std::cout << m << std::endl ;
```

The same can be achieved by using the standalone application SequenceMatrixCreator which is a wrapper around an instance of SequenceMatrixCreator.

```
SequenceMatrixCreator --bed /some/bed/file.bed --fasta /some/fasta/file.fasta --from -100 --to 100
```

For help, enter :

```
SequenceMatrixCreator -h
```


# 3 Partitioning algorithms

Once we have matrices containing read densities, or DNA sequences, from regions of interest, we may be interested in partitioning them in order to isolate leading trends either in terms of chromatin features (for read densities) or in terms of motif contents (for DNA sequences).

Each partitioning algorithm is relies on a common expectaction-maximization (EM) core that was described in ChIPPartitioning algorithm (here refered to as EMRead) in https://academic.oup.com/bioinformatics/article/30/17/2406/2748187.

Because all of these algorithms are derived from ChIPPartitioning/EMRead, it is strongly advised to read to corresponding section. The details presented in that section can be necessary to understand the other algorithms.

![](images/em.png?raw=true)

**EMRead and EMSequence algorithms**
On the right, EMRead algorithm. From top to bottom : the input is a 2D matrix of read densities, with one region per row. The regions are divided in S slices of lengh L'=L-S+1 where L is the matrix number of columns (bins) and S is the shifting freedom. Each slice begin at offset 1, 2, ..., S. The K different classes are represented by a vector of size L' containing the weighted aggregation of signal in each class. The partition is optimized using an expectation-maximization algorithm. During the E-step, for each slice, a likelihood to belong to each class is computed. This is equivalent to scanning the entire regions of length L with the models of length L'. Given the class probability, which can be interpreted as the class size, a posterior probability is computed. The posterior probabilities are interpreted as class membership scores. During the M-step, the class models are updated using the posterior probabilities as weights for the aggregation procedure. The class probabilities are also updated. The E- and M-steps are then repeated iteratively for a given number of times.
On the left, EMSequence algorithm. From top to bottom : the input is a 2D matrix of DNA sequences encoded as integers usin A->0, C->1, G->2, T->3, N->4. The sequences are divided in S slices of lengh L'=L-S+1 where L is the matrix number of columns (bins) and S is the shifting freedom. Each slice begin at offset 1, 2, ..., S. The K different class are represented by a vector of size L' containing the weighted aggregation of signal in each class. The K different classes are represented by a letter probability matrix (LPM) of dimension L'x4 (for A, C, G, T) that contain the probability of each base at each position in the model. Then the E- and M-steps are strictly similar to EMRead. At the end of the M-step, the class models are also updated as a weighted aggregation. However, this results in creation of a 2D matrix (a LPM) instead of a vector.


## 3.1 EMRead (ChIPPartitioning)

EMRead is a probabilistic partitioning algorithm that softly clusters a set of genomic regions based on their signal shape (as opposed to the absolute values) resemblance. To ensure proper comparisons between the regions, the algorithm allows to offset one region compare to the other to retrieve a similar signal at different offsets and to flip the signal orientation. Finally, it has been demonstrated to be really robust to sparse data. 

EMRead (a graphic representation of the algorithm can be found in the figure above, in the left panel) models the read density signal over N region of length L as having being sampled from a mixture of K different read density models (classes), using L independent Poisson distributions for each position. The number of reads sequenced over this region is then the result of this sampling process. Each class model is represented by a vector of size L'=L-S+1 where S is the number of shift states allowed. 

For now, let us assume S=1. In this case, L'=L and the class models are vectors as long as the regions. Because EMRead is a probabilistic partitioning program, it assigns all regions to all classes, with different probabilities (that can be interpreted as weights). A class model is then simply the mean signal present over all regions, weighted by the probability that each region belong to that class. Thus, a class model is a weighted aggregation.

Now, setting S>1 allows to realign the regions in order to retrieve a similar signal, in two regions, but at different position in each regions. For instance, the 1st region may have a strong signal at its beginning and the 2nd may have a strong signal at its end. A regular position-wise comparison of the signal will fail to detect that both regions have a strong signal stretch and will even lead to the conclusion that both regions are dissimilar. More formally, this is done by subdividing each regions in slices of length L'=L-S+1 that starts at each possible offset 1, 2,...,S. Then, instead of assigning the entire regions (rows) to classes, EMRead assigned the slices to classes. The class models are then simply the aggregation of the slices, weighted by the probability that these slices were assigned to the classes. 

The alignment process can further be improved by allowing flip. In this case, an addition flip state (reverse) are allowed, instead of forward only. In this case, each region (if S=1) or slice (if S>1), are compared with the class models in forward and reverse orientation. For the 2nd case, the region/slice is simply flipped. That is, it is from the end to the start )the 1st position becomes the last and so one). The class model computation procedure stays the same, except that the reverse regions/slices are also taken into account during the computation of the weighted aggregations.

EMRead is an iterative (it is an EM) algorithm. The class models are initialized by randomly assigning the regions to the classes. This also initialize the overall class probability (which can be interpreted as the sizes/importances of the classes). During the M-step, each region/slice is compared with each class model and likelihood are computed. Given the class probabilities, the posterior probabilities corresponding to the posterior probabilities are computed. Once the posterior probabilities are known, the model classes are updated. Eventually, this process is repeated for a given number of iterations and the data partition defined by the posterior probabilities.

### 3.1.1 Background class

A background class can be included to model noise. A background class has the following properties :


1. Its model has the same length as all the other classes models (L').
2. Its model is initialised with the mean number of reads in the entire dataset at every position (it is a flat signal).
3. Its model is never updated, it remains stable during the entire partitioning process.
4. It is always the last class. In a 6 class partition, classes 1 to 5 are regular classes and the 6th is the background class.

Except these specificities, the computations remains the same. Each region is computed a likelihood given this model and is attributed a probability of belonging to this class at each iteration.

### 3.1.2 Example

Given a read density matrix, a partition can be performed using the EMRead class in Clustering. The following snippet partitions a dataset into 2 to 10 classes. Each partition is optimized for 20 iterations and a shifting of 21 bins is allowed together with flipping. The computations will be parallelized over 4 threads.

```
#include <iostream>
#include <vector>
#include "Random/Random.hpp"   // rand_string()
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Clustering/EMRead.hpp"

// data is Matrix2D<int>

// will contain all the partitions 
vector<Matrix4D<double>> partitions(9) ;

// generates a random seed of 30 characters long
// to seed the random number generators
std::string seed = rand_string(30) ;

size_t n_iter     = 20 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

for(size_t k=2; k<11; i++)
{
	EMRead em(data,
        	  k,
        	  n_iter,
        	  n_shift,
        	  flip,
        	  bckg_class,
        	  seed,
        	  n_threads) ;
	em.classify() ;
	// keeps the posterior probabilities in memory
	partitions[k-2] = em.get_post_prob()
	// dump the posterior probabilities in a binary 
	// file named with the seed
	partitions[k-2].save(seed) ;
	// display the dimensions of the probability matrix
	// should be equal to data.get_nrow(), k, n_shift, flip+1
	std::cout << partitions[k-2].get_dim() << std::endl ;
}
```
The output is a 4D matrix of doubles with dimensions NxKxSxF where N is the number of regions, K the number of classes, S the number of shift states and F the number of flip states (1 normally, 2 with flip). This matrix contains the probability that each slice (3rd dimension) of each region (1st dimensions) in a given orientation (4th dimension) belongs to each class (2nd dimension). For the 4th dimensions, the index 0 always relate to the "forward" orientation and the index 1 to the "reverse" orientation. If no flip is allowed, the 4th dimension has length 1, with flip it has length 2.

The same can be achieved by using the standalone application EMRead which is a wrapper around an instance of EMRead. In that case, the 4D probability matrix is dumped in a file in binary format. 

```
EMRead --read /path/to/data.mat --class k --shift 21 --flip --iter 20 --seed tralalala --thread 4 --out /path/to/prob.mat4d
```
For help enter :

```
EMRead -h
```

Loading a 4D matrix can be done using :

```
#include "Matrix/Matrix4D.hpp"

Matrix4D<double> prob ;
prob.load("/path/to/binary_file") ;

```

### 3.1.3 Extracting the class models

The above examples showed how to access the posterior probabilities rather than the computed models. Accessing the models from an EMRead instance is possible once the classify() method has been called. The models are stored in a 3D matrix of dimensions KxL'x1 where K is the number of classes and L' is the length of the models (L'=L-S+1).

```
#include <iostream>
#include "Clustering/EMRead.hpp"
#include "Matrix/Matrix3D.hpp"
EMRead em(...) ;
em.classify() ;

// extract the models
Matrix3D<double> models = em.get_read_models() ;
// display the models
std::cout << models << std::endl ;

```

## 3.2 EMSequence

EMSequence is a probabilistic partitioning algorithm that softly clusters a set of genomic regions based on their DNA sequences motif content (a graphic representation of the algorithm can be found in the figure above, in the right panel). EMSequence algorithm has been derived from EMRead algorithm to handle DNA sequences. In essence, EMSequence is a de novo motif finding algorithm. All the computations are strictly the same with the exceptions of everything related to the class models. The differences are :

1. the class models are not vectors but 2D matrices that contain letter probability matrices (LPMs)
2. the data likelihood computation do not depend on Poisson distribution anymore but on the base distributions described in the LPMs
3. the class model update procedure has been modified and the weighted aggregation of DNA sequences result in the generation of LPMs.

Let us cover the above mentioned points in more details. 

First, the class models are LPMs. An LPM is a matrix of dimensions L'x4 that describe a consensus sequence. An LPM contains the probability of A, C, G and T respectively at each of the L' positions modeled. For instance, the following LPM

| A    | C    | G    | T    |
| ---- |----- | ---- | -----|
| 1.00 | 0.00 | 0.00 | 0.00 |
| 0.00 | 0.50 | 0.50 | 0.00 |
| 0.33 | 0.00 | 0.33 | 0.33 |

describes a DNA sequence of length 3 that bears always bear A in 1st position, C or G (equally likely) in 2nd position and A or G or T (equally likely) in 3rd position. The rows always sum to 1 as the values are per-position probabilies. Also, LPM models are 0-order (or mono-nucleotide) model. The occurence of a base at a given position has no impact any on the probability of occurence of any other nucleotide at any other position. Also, bear in mind that 0 values are not authorized in LPMs for two reasons : i) statistically, a probability of 0 is probably due to a sampling error and ii) computationally, to avoid divisions by 0. Thus, a minimal probability value is enforced in LPMs and in any probability computation.

Second, the likelihood of a given sequence of length L' given a model of dimension L'x4 is computed as the product of the probability of appearance of each base at each position. Each row of the model thus contains a discrete probability distribution.

Third, compared to EMRead,the model class update has been modified. In EMRead, the model update performs a weighted aggregation of the signal, for each class, using the assignment probabilies as weights. Obviously, DNA characters cannot be aggregated the same. An expression such as "0.5A + 0.25C", instead of giving a single value, results in the following probability distribution "0.5/(0.5+0.25)A, 0.25/(0.5+0.25)C, 0.0/(0.5+0.25)G, 0.0/(0.5+0.25)T" which is "0.66A, 0.33C, 0.0G, 0.0T". Eventually, aggegating all the sequence slices, weighted by their probabilities, results in the creation of one LPM per class.

Otherwise, the overall algorithm design is the same as for EMRead. The partition is iteratively optimized using an EM scheme. The final output is also strictly identical to EMRead. It is a 4D matrix of doubles with dimensions NxKxSxF where N is the number of sequences, K the number of classes, S the number of shift states and F the number of flip states (1 normally, 2 with flip). This matrix contains the probability that each slice (3rd dimension) of each sequence (1st dimensions) in a given orientation (4th dimension) belongs to each class (2nd dimension). For the 4th dimensions, the index 0 always relate to the "forward" orientation and the index 1 to the "reverse" orientation. If no flip is allowed, the 4th dimension has length 1, with flip it has length 2.

### 3.2.1 Background class

A background class can be included, exactly as in EMRead. The background class model is initialised with the average base probabilities (over the entire dataset) at every position

### 3.2.2 Example

The following snippet of code uses a EMSequence instance, from Clustering, to partition a dataset of DNA sequences into 6 classes, for 200 iterations, allowing a shifting freedom of 21 and flipping. The computations will be parallelized over 4 threads. Note that the input has to be an integer matrix containing one nucleotide per cell using the following encoding : A->0, C->1, G->2, T->3, N->4. If the input is a character matrix, a conversion can be performed using the dna::char_to_int() function from dna_utility.cpp in Utility.

```
#include <string>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Utility/dna_utility.hpp"   // for char_to_int()
#include "Random/Random.hpp"         // rand_string()
#include "Clustering/EMSequence.hpp"

// generates a random seed of 30 characters long
// to seed the random number generators
std::string seed = rand_string(30) ;

// ====================== case with int matrix ======================

// data_int is a Matrix2D<int>, no conversion needed

size_t n_class    = 6 ;
size_t n_iter     = 200 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

EMSequence em(data_int,
              n_class,
              n_iter,
              n_shift,
              flip,
              bckg_class,
              seed,
              n_threads) ;
em.classify() ;
// the posterior probabilities in memory
Matrix4D<double> prob = em.get_post_prob()
// dump the posterior probabilities in a binary 
// file named with the seed
prob.save(seed) ;

// ====================== case with char matrix ======================

// data_char is a Matrix2D<char>, conversion needed

em = EMSequence(dna::char_to_int(data_char),
                n_class,
                n_iter,
                n_shift,
                flip,
                bckg_class,
                seed,
                n_threads) ;
em.classify() ;
// the posterior probabilities
Matrix4D<double> prob2 = em.get_post_prob()
// dump the posterior probabilities in a binary 
// file named with the seed
prob2.save(seed) ;
```

The same can be achieved by using the standalone application EMSequence which is a wrapper around an instance of EMSequence. 

```
EMSequence --seq /path/to/data.fasta --class 6 --shift 21 --flip --iter 200 --seed poupoum --thread 4 --out /path/to/prob.mat4d
```

For help enter :

```
EMSequence -h
```

As for EMRead, the 4D probability matrix is dumped in a file in binary format. Loading the matrix can be done using :

```
#include "Matrix/Matrix4D.hpp"

Matrix4D<double> prob ;
prob.load("/path/to/binary_file") ;
```

### 3.2.3 Extracting the class models

The above examples showed how to access the posterior probabilities rather than the computed models. Accessing the models from a EMSequence instance is possible once the classify() method has been called. The models are stored in a 3D matrix of dimensions KxL'x4 where K is the number of classes and L' is the length of the models (L'=L-S+1) and 4 for A, C, G and T respectively.

```
#include <iostream>
#include "Clustering/EMSequence.hpp"
#include "Matrix/Matrix3D.hpp"
EMSequence em(...) ;
em.classify() ;

// extract the models
Matrix3D<double> models = em.get_sequence_models() ;
// display the models
std::cout << models << std::endl ;
```

## 3.3 EMConsensusSequence

Imagine for a minute that instead of having regular DNA sequences (for instance 'A,C,G'), we have consensus sequences (the rational will become clearer lower). That is, at each position, several different nucleotides can be present (for instance 'A,C/G,A/G/T). Consensus sequences are represented using LPMs. A consensus sequence of length L is represented using a LPM of dimensions Lx4. The above case would be represented as follow :


| A    | C    | G    | T    |
| ---- |----- | ---- | -----|
| 1.00 | 0.00 | 0.00 | 0.00 |
| 0.00 | 0.50 | 0.50 | 0.00 |
| 0.33 | 0.00 | 0.33 | 0.33 |

This is exactly what happens when EMSequence class models are extracted. After a partitioning of N DNA sequences with EMSequence into K classes, K LPMs can be extracted. If for some reason partitioning these LPMs becomes of interest, EMConsensusSequence should be used. The algorithm works exactly as EMSequence except that it works on LPMs (if you think of it, a DNA sequence is a special case of LPM).

The only difference with EMSequence is that, instead of taking a 2D matrix as input, EMConsensusSequence takes a 3D matrix as input. Its dimensions are NxLx4 where N is the number of consensus sequences, L is the length of the consensus sequences and 4 for A, C, G and T respectively.

The class models are LPMs as for EMSequence, as aggregating LPMs is equivalent to computing a weighted aggregation of DNA sequences, as in EMSequence.

The output format is identical to EMRead and EMSequence.


### 3.3.1 Background class

A background class can be included, exactly as in EMSequence. The background class model is initialised with the average base probabilities (over the entire dataset) at every position

### 3.3.2 Example

The following snippet of code uses a EMConsensusSequence instance, from Clustering, to partition a set of DNA consensus sequences into 4 classes, for 200 iterations, allowing a shifting freedom of 21 and flipping. The computations will be parallelized over 4 threads.

```
#include <string>
#include "Matrix/Matrix3D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Random/Random.hpp"                   // rand_string()
#include "Clustering/EMConsensusSequence.hpp"

// generates a random seed of 30 characters long
// to seed the random number generators
std::string seed = rand_string(30) ;

// data_seq is a Matrix3D<double> containing
// the consensus sequences


size_t n_class    = 4 ;
size_t n_iter     = 200 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

EMConsensusSequence em(data_seq,
                       n_class,
                       n_iter,
                       n_shift,
                       flip,
                       bckg_class,
                       seed,
                       n_threads) ;
em.classify() ;
// the posterior probabilities in memory
Matrix4D<double> prob = em.get_post_prob()
// dump the posterior probabilities in a binary 
// file named with the seed
prob.save(seed) ;
```

The same can be achieved by using the standalone application EMConsensusSequence which is a wrapper around an instance of EMConsensusSequence. 

```
EMConsensusSequence --consseq /path/to/data.mat3d --class 4 --shift 21 --flip --iter 200 --seed tsointsoin --thread 4 --out /path/to/prob.mat4d
```

For help enter :

```
EMConsensusSequence -h
```

As for EMRead and EMSequence the 4D probability matrix is dumped in a file in binary format. Loading the matrix can be done using :

```
#include "Matrix/Matrix4D.hpp"

Matrix4D<double> prob ;
prob.load("/path/to/binary_file") ;
```

### 3.3.3 Extracting the class models

For EMConsensusSequence instances, the models are stored in a 3D matrix of dimensions KxL'x4, as for the EMSequence class.

```
#include <iostream>
#include "Clustering/EMConsensusSequence.hpp"
#include "Matrix/Matrix3D.hpp"
EMConsensusSequence em(...) ;
em.classify() ;

// extract the models
Matrix3D<double> models = em.get_sequence_models() ;
// display the models
std::cout << models << std::endl ;
```


## 3.4 EMJoint

EMJoint is the generalization of EMRead, EMSequence and EMConsensusSequence. It can perform i) the partition of a set of regions based on on or more read densities for the same genomic regions or ii) one or more read densities for the same genomic regions and their DNA sequences. The DNA sequence can be regular DNA sequences (stored in a 2D matrix) or consensus sequences (stored in a 3D matrix).


EMJoint takes as input several matrices (called data layers) that describe different signal over strictly identical genomic regions. For instance, if a read density matrix R and a sequence matrix S are given, then S(100,382) contains the base that is found in the sequence of the 101th regions, at position 383. Consequently, R(100,382) must contain the number of reads found at the 383th position of the 101th region. This also implies that the bin size of the different matrices should be the same. If a sequence matrix is given, then all matrix bin sizes must be 1bp (each cell in a sequence matrix contains a 1bp sequence). If only read density matrices are given, the 1bp rule does not apply but the bin size should be the same among the matrices.

EMJoint uses jointly the framework of EMRead for read density related computations and of EMSequence or EMConsensusSequence for sequence related computations. The only things that differ from EMRead or EMSequence/EMConsensusSequence are 

1. the likelihood computations
2. how the model classes are represented

Regarding the likelihood computations, a joint likelihood is computed that encompasses all the data layers. Let us assume that a read density and a sequence layer have been given. The likelihood of a region is computed as the product of the individual layer likelihood (using there dedicated models). For read density signal, the computations are defined in EMRead and for the read layer(s), in EMSequence for the sequence layer and in EMConsensusSequence for consensus sequences. To be precise, these computations are implemented in the ReadLayer, SequenceLayer and ConsensusSequenceLayer classes in Clustering (these classes implements the data storage, the class models and the likelihood computations).

Regarding the model class representations each class is represented by several distinct signal models. There are exactly 1 model per class per signal layer. Read density layers are modeled using signal vectors, as in EMRead, that are the weighted aggregation of the read densities in a class. Sequence layers are modeled using LPMs, as in EMSequence. Let us take the example of a case with 2 read layers and 1 sequence layer and 3 classes. In this case, there will be 3*2=6 read density models and 3*1 LPMs. Each class will contain 2 read density models and 1 LPM.

Otherwise, the computations are strictly identical to what EMRead and EMSequence do, including the output format.

### 3.4.1 Example 1

The following snippet of code uses a EMJoint instance, from Clustering, to partition a set of genomic regions for which there are read density data and their DNA sequences, into 6 classes, for 20 iterations, allowing a shifting freedom of 21 and flipping. The computations will be parallelized over 4 threads. Note that the sequence matrix has to be an integer matrix containing one nucleotide per cell using the following encoding : A->0, C->1, G->2, T->3, N->4.


```
#include <cassert>                   // assert(() 
#include <vector>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Random/Random.hpp"         // rand_string()
#include "Clustering/EMJoint.hpp"

// generates a random seed of 30 characters long
// to seed the random number generators
std::string seed = rand_string(30) ;

// data_seq is a Matrix2D<int> containing
// sequences in int format

// data_read is a Matrix2D<int> containing
// read densities

// matrix dimensions should perfectly fit
assert(data_seq.get_nrow() == data_read.get_nrow() and
       data_seq.get_ncol() == data_read.get_ncol()) ;

size_t n_class    = 6 ;
size_t n_iter     = 20 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

std::vector<Matrix2D<int>> vect_read ;
vect_read.push_back(data_read) ;

EMJoint em(vect_read,
           data_seq,
           n_class,
           n_iter,
           n_shift,
           flip,
           bckg_class,
           seed,
           n_threads) ;
em.classify() ;
// the posterior probabilities
Matrix4D<double> prob = em.get_post_prob() ;
// dump the posterior probabilities in a binary 
// file named with the seed
prob2.save(seed) ;
```
The same can be achieved by using the standalone application EMJoint which is a wrapper around an instance of EMJoint. 

```
EMJoint --read /path/to/dnase.mat --seq /path/to/seq.mat --class 6 --shift 21 --flip --iter 20 --seed guillaume --thread 4 --out /path/to/prob.mat4d
```

For help, enter :

```
EMJoint -h
```

As for EMRead and EMSequence, the 4D probability matrix is dumped in a file in binary format. Loading the matrix can be done using :

```
#include "Matrix/Matrix4D.hpp"

Matrix4D<double> prob ;
prob.load("/path/to/binary_file") ;
```

### 3.4.2 Example 2

The following snippet of code uses a EMJoint instance, from Clustering, to partition a set of genomic regions for which there are twp read density datasets - let us say from a DNase-seq and a ChIP-seq experiments - into 6 classes, for 20 iterations, allowing a shifting freedom of 21 and flipping. The computations will be parallelized over 4 threads.

```
#include <cassert>                   // assert()
#include <vector>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Random/Random.hpp"         // rand_string()
#include "Clustering/EMJoint.hpp"

// generates a random seed of 30 characters long
// to seed the random number generators
std::string seed = rand_string(30) ;

// dnase_read is a Matrix2D<int> containing
// DNase-seq read densities

// chipseq_read is a Matrix2D<int> containing
// ChIP-seq read densities

// matrix dimensions should perfectly fit
assert(dnase_read.get_nrow() == chipseq_read.get_nrow() and
       dnase_read.get_ncol() == chipseq_read.get_ncol()) ;

size_t n_class    = 6 ;
size_t n_iter     = 20 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

std::vector<Matrix2D<int>> vect_read ;
vect_read.push_back(dnase_read) ;
vect_read.push_back(chipseq_read) ;

EMJoint em(vect_read,
           n_class,
           n_iter,
           n_shift,
           flip,
           bckg_class,
           seed,
           n_threads) ;
em.classify() ;
// the posterior probabilities
Matrix4D<double> prob = em.get_post_prob() ;
// dump the posterior probabilities in a binary 
// file named with the seed
prob2.save(seed) ;
```
The same can be achieved by using the standalone application EMJoint.

```
EMJoint --read /path/to/dnase.mat,/path/to/chipseq.mat --class 6 --shift 21 --flip --iter 20 --seed meurice --thread 4 --out /path/to/prob.mat4d
```

### 3.4.3 Extracting the class models

For EMJoint instances, the read models are stored in a 3D matrix of dimensions KxL'x1 as for EMRead. The sequence models are stored in 3D matrix of dimensions KxL'x4 as for EMSequence. However, two different methods have to be called depending on the DNA sequence types. If the sequence are regulare DNA sequences, then EMJoint::get_sequence_models() should be used. If consensus sequences were given, EMJoint::get_consensus_sequence_models() should be used.

```
#include <iostream>
#include "Clustering/EMJoint.hpp"
#include "Matrix/Matrix3D.hpp"
EMConsensusSequence em(...) ;
em.classify() ;

// ============= case with regular DNA sequences =============

// extract the models
Matrix3D<double> models = em.get_sequence_models() ;
// display the models
std::cout << models << std::endl ;

// ============= case with consensus sequences =============

// extract the models
Matrix3D<double> models = em.get_consensus_sequence_models() ;
// display the models
std::cout << models << std::endl ;

```

## 3.5 Remark about multithreading

All the above examples have distributed the computations over different threads. Multithreading is not always desirable. These classes implementing the partitioning algorithms implement an "embarassingly parallel" distribution scheme. The total workload is equally distributed over the different threads. However, if your system does not have enought free CPUs or if the amount of computations per thread (roughly number of regions /  number of threads) is too low, multithreading may not be beneficial.


# 4 Post processing operations

With a given data partition, that is with a matrix of posterior probabilities, several desirable operations can be performed such :

1. Realigning a dataset given a partition. This partition may have been performed using this same dataset, in which case is equivalent as extracting the class models. If the partition has been performed on another dataset, this allows to co-visualize different layers of data.
2. Extending a set of models of length L' to a new length L''=L'+E with the original model as centers of the extended ones.
3. Extracting the data assigned to a class.

## 4.1 Realigning a dataset given a partition

![](images/realignment.png?raw=true)
**Data realignment process**
Toy example showing how to realign data to compute class models. 
All the partitioning algorithms use the same overall procedure to compute/update their class models. In brief, all the slices are aggregated together, using the posterior probabilities as weights. This procedure is somewhat equivalent to a gap-free realignment procedure. Thus, given a set of posterior probabilities, it is possible to realign another dataset to compute the corresponding aggregations.
Toy example illustrating the data realignment procedure. All partitioning algorithms use the same overall procedure to compute/update their class models. In brief, all the slices are aggregated together, using the posterior probabilities as weights. This procedure is somewhat equivalent to a gap-free realignment procedure. In the center, a dataset containing DNA sequences have been partitioned using EMSequence. One of the classes - the one displayed - has identified the 'AGGC' subsequence that is present in all the sequences. The slice alignment and the corresponding class aggregation are displayed. 
On the right, a data read density matrix, representing strictly the same regions as the sequence matrix, is realigned using the posterior probabilities computed by EMSequence. The slice alignment and their agrgegation is displayed. Co-vizualizing both models allow to visualize a read density over a sequence motif. Here the motif, displayed as a DNA logo, is highly conserved and presents a high read coverage.
On the edges, the C++ classes that implemente these procedures, their inputs and their outputs are displayed for read density and DNA sequence data. These C++ classes compute always compute the models of all classes at once an return them as 2D matrix that contains both the models (read density vectors or LPMs, indicated in purple) and the class probabilities (indicated in grey-blue). 

This procedure can be used in two different cases.

1. To recompute class models. It is possible to extract the class models directly from the instances that implement a partitioning proceudre, such as an EMRead instance using EMRead::get_read_models(). Additionally, tt is also possible to recompute the class models given the data and the posterior probability.

2. To align a dataset B given a partition made on a dataset A. It is absolutely feasible to run a partitioning on a given matrix A, for instance DNA sequences, using EMSequence, and to subsequently use the obtained posterior probabilities to compute the class models, using another data matrix, let us say B of DNase-seq reads. This procedure allows to realign a dataset B as A in order to visualize how different types of signals co-occure. The only things that should be taken care of is that matrices A and B should have the same dimensions and that the genomic positions inside both matrices are strictly identical. The toy example displayed in the figure above depicts a case in which a set of regions, for which there are the DNA sequences and a DNase-seq indicating chromatin accessibility. Partitioning the regions based on their DNA sequences allows to find 'AGGC' as prefered the binding sequence for the TF. Additionally, the DNase-seq read densities over these regions are available. Vizualizing the read density centered on the TF motif may be interesting. Eventually, realigning the read densities based on the DNA sequence alignment allows to vizualize them together by overlapping the read density and the sequence motifs. The motif has a constant read density meaning it seems somewhat accessible.


### 4.1.1 Realigning read densities

This can be done using an instance of the ReadModelComputer class from Clustering. The output is a 2D matrix with dimensions KxL'+1 where K is the number of classes and L' the model length. The 1st column always contains the class probabilities. The following columns contains the class models, as depicted in the figure above.

The following code snippet shows how to recompute the read class models from a partition (case 1 above).

```
#include <iostream>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix3D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Clustering/EMRead.hpp"
#include "Clustering/ReadModelComputer.hpp"

// data is Matrix2D<int> containing read densities

bool bckg_class  = false ;
size_t n_threads = 4 ;

// create the partition
EMRead em(data,
          ...,
          bckg_class,
          n_threads) ;
em.classify() ;

// extract the probabilities
Matrix4D<double> prob = em.get_post_prob()

// recompute the models
ReadModelComputer rmc(data,
                      prob,
                      bckg_class,
                      n_threads) ;
Matrix2D<double> model = rmc.get_model() ;

// display class probabilities
for(size_t i=0; i<model.get_nrow(); i++)
{	std::cout << "prob class " << i+1 << " " << model(i,0) << std::endl ; }

// display class models
for(size_t i=0; i<model.get_nrow(); i++)
{	std::cout << "model class " << i+1 << std::endl ;	
	for(size_t j=1; j<model.get_ncol(); j++)
	{	std::cout << model(i,j) << " " ; }
	std::cout << std::endl ;
}

```

The same can be done using the autonomous application ProbToModel :

```
ProbToModel --read /path/to/read.mat -- prob /path/to/prob.mat4d --thread 4
```

For help, enter :

```
ProbToModel -h
```

### 4.1.2 Realigning sequences

This can be done using an instance of the SequenceModelComputer class from Clustering. The output is a 2D matrix with dimensions (4xK)x(L'+1) where K is the number of classes and L' the model length. The 1st column always contains the class probabilities. The following columns contains the individual class LPMs, as depicted in the figure above. Each class LPM is stored on 4 rows. For instance, with 3 classes, rows 1-4 contains class 1, rows 5-8 class 2 and rows 9-12 class 3.

The following code snippet show how to realign DNA sequences given a partition done on another dataset, for instance on a read density dataset (case 2 above).

```
#include <iostream>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix3D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Clustering/EMRead.hpp"
#include "Clustering/SequenceModelComputer.hpp"

// data_read is Matrix2D<int> containing read densities

// data_seq  is Matrix2D<int> containing DNA sequences
// encoded as A->0, C->1, G->2, T->3

bool bckg_class  = false ;
size_t n_threads = 4 ;

// create the partition
EMRead em(data_read,
          ...,
          bckg_class,
          n_threads) ;
em.classify() ;

// extract the probabilities
Matrix4D<double> prob_read = em.get_post_prob()

// realign the DNA sequences given the
// read density partition
SequenceModelComputer smc(data_seq,
                          prob_read,
                          bckg_class,
                          n_threads) ;
Matrix2D<double> model = smc.get_model() ;

// display class probabilities
size_t n_class = model.get_nrow() / 4 ;
for(size_t i=0; i<n_class; i++)
{	std::cout << "prob class " << i+1 << " " << model(i*4,0) << std::endl ; }

// display class models
for(size_t i=0; i<n_class; i++)
{	std::cout << "model class " << i+1 << std::endl ;	
	for(size_t i_2=0; i_2<4; i_2++)
	{	for(size_t j=1; j<model.get_ncol(); j++)
		{	std::cout << model(i+i_2,j) << " " ; }
		std::cout << std::endl ;
	}
	std::cout << std::endl ;
}

```

The same can be done using the autonomous application ProbToModel :

```
ProbToModel --seq /path/to/seq.mat --prob /path/to/prob.mat4d --thread 4
```

### 4.1.3 Realigning consensus sequences

This can be done using an instance of the ConsensusSequenceModelComputer class from Clustering. The output is, as for SequenceModelComputer, a 2D matrix with dimensions (4xK)x(L'+1) where K is the number of classes and L' the model length. The 1st column always contains the class probabilities. The following columns contains the individual class LPMs, as depicted in the figure above. Each class LPM is stored on 4 rows. For instance, with 3 classes, rows 1-4 contains class 1, rows 5-8 class 2 and rows 9-12 class 3.

The following code snippet show how to realign consensus sequences given a partition done on another dataset, for instance on a read density dataset (case 2 above).

```
#include <iostream>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix3D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "Clustering/EMRead.hpp"
#include "Clustering/ConsensusSequenceModelComputer.hpp"

// data_read is Matrix2D<int> containing read densities

// data_seq  is Matrix3D<double> containing consensus sequences
// encoded as A->0, C->1, G->2, T->3

bool bckg_class  = false ;
size_t n_threads = 4 ;

// create the partition
EMRead em(data_read,
          ...,
          bckg_class,
          n_threads) ;
em.classify() ;

// extract the probabilities
Matrix4D<double> prob_read = em.get_post_prob()

// realign the DNA sequences given the
// read density partition
ConsensusSequenceModelComputer smc(data_seq,
                                   prob_read,
                                   bckg_class,
                                   n_threads) ;
Matrix2D<double> model = smc.get_model() ;
```

The same can be done using the autonomous application ProbToModel :

```
ProbToModel --consseq /path/to/consseq.mat3d --prob /path/to/prob.mat4d --thread 4
```


## 4.2 Extending class models

![](images/extension.png?raw=true)
**Model extension procedure**
Toy example showing the model extansion procedure by extending the models obtained in the previous figure.
The models need to be extended by E=4 columns overall. That is, E/2=2 columns will be added on each side. For this, the original data matrix (let us call it the "short matrix") from which the original class models were computed (let us call them the "short models") is extended. The data realignment is performed as shown in the previous figure, with the posterior probabilities that have been used to compute the short model. Because the data matrix has been widen, the slices are longer, which results in the creation of longer class models.
On the edges, the C++ classes that implemente these procedures, their inputs and their outputs are displayed for read density and DNA sequence data. These C++ classes compute always compute the models of all classes at once an return them as 2D matrix that contains both the models (read density vectors or LPMs, indicated in purple) and the class probabilities (indicated in grey-blue). 

The above figure allowed to co-visualize two models together. The sequence model shows a core motif of 4bp (AGGC) that is ultra conserved but does not provide any information about the flanking regions. Additionally, the motif indicates an opening through the read density aggregation. However, we don't know how open is the motif compared to the flanking regions. Is it more or less open? And what about the motif? Is it a short conserved motif to which a transcription factor (TF) binds or is it a short stretch inside a longer conserved region, such as a repeated element?

In order to visualize the flanking regions, the class models can be extended, based on the data.

In this case, we have partitioned a DNA sequence matrix S of dimensions NxL using K classes, with a shifting freedom S and with flipping. The posterior matrix probability P has dimensions NxKxSx2 (region, class, shift, flip) and the K models have a length of L'=L-S+1.

Extending the models is obtained by computing a larger matrix S<sup>ext</sup> of dimensions NxL'' where L''=L+E. In this case E is the extra number of columns to add. Care should be taken to construct S<sup>ext</sup> by adding exactly E/2 columns one each side of S such that S is contained in the central part of S<sup>ext</sup>.

Once S<sup>ext</sup> has been constructed, the class model computation step (here for DNA sequence data) can be applied on it using the posterior probability matrix P. This will results in the creation of K models of length L''-S+1. The K original models will be contained in the central part of the K extended models. Extending read models is done exactly the same, excepted that the class model computation step for read densities should be used.

In the case aboce, the read and sequence models were extended by E=4 columns (2 on each side). This revealed two things. First, the sequence motif has an extra conserved A. Second, the motif accessibility in in fact weaker than its direct flanking regions which indicates a protection. This pattern is compatible with a footprint, suggesting that this is a class representing an occupied TF binding site.

### 4.2.1 Extending read density models

The following snippet of code indicates how to extend read density models.

```
#include <iostream>
#include <string>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "GenomicTools/CorrelationMatrixCreator.hpp"
#include "Clustering/EMRead.hpp"
#include "Clustering/ReadModelComputer.hpp"

std::string bed_path = "/path/to/file.bed" ;
std::string bam_path = "/path/to/file.bam" ;
std::string bai_path = "/path/to/file.bai" ;

// original matrix creation and partitioning
int from     = -1000 ;
int to       =  1000 ;
int bin_size =     1 ;
std::string method = "read" ;

CorrelationMatrixCreator mc(bed_path,
                            bam_path,
                            bai_path,
                            from,
                            to,
                            bin_size,
                            method) ;
Matrix2D<int> data = mc.create_matrix() ;

size_t n_class    = 4 ;
size_t n_iter     = 20 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

EMRead em(data,
          n_class,
          n_iter,
          n_shift,
          flip,
          bckg_class,
          std::string("the_bird_is_a_word"),
          n_threads) ;
em.classify() ;
Matrix4D<double> prob = em.get_post_prob() ;
Matrix3D<double> model = em.get_read_models() ;


// procedure to extend the models
// 1 modify from/to range
int   n_ext       = 1000*bin_size ; // add 500 columns of 1bp on each side
int   n_ext_right = n_ext/2 ;
int   n_ext_left  = n_ext - n_ext_right ;
from -= ext_left ;
to   += ext_right ;

// 2 create the extended matrix
CorrelationMatrixCreator mc2(bed_path,
                             bam_path,
                             bai_path,
                             from,
                             to,
                             bin_size,
                             method) ;
Matrix2D<int> data_ext = mc2.create_matrix() ;

// 3 compute the models
ReadModelComputer rmc(data_ext,
                      prob,
                      bckg_class,
                      n_threads) ;
Matrix2D<double> model_ext = rmc.get_model() ;

// 4 display the extended models
std::cout << model_ext << std::endl ;

```

The same could be performed using the autonomous application ReadModelExtender that implements this procedure.

```
ReadModelExtender --bed /path/to/file.bed --bam /path/to/file.bam --bai /path/to/file.bai --prob /path/to/prob.mat4d --from -1000 --to 1000 --binSize 1 --ext 1000 --method read --thread 4
```

For help, enter :

```
ReadModelExtender -h
```


### 4.2.2 Extending DNA sequence models

The following snippet of code indicates how to extend DNA sequence models.

```
#include <iostream>
#include <string>
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"
#include "GenomicTools/SequenceMatrixCreator.hpp"
#include "Clustering/EMSequence.hpp"
#include "Clustering/SequenceModelComputer.hpp"

std::string seq_path = "/path/to/file.fasta" ;
std::string bed_path = "/path/to/file.bed" ;

// original matrix creation and partitioning
int from     = -1000 ;
int to       =  1000 ;

SequenceMatrixCreator mc(bed_path,
                         seq_path,
                         from,
                         to) ;
Matrix2D<int> data = mc.create_matrix() ;

size_t n_class    = 4 ;
size_t n_iter     = 20 ;
size_t n_shift    = 971 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;

EMSequence em(data,
              n_class,
              n_iter,
              n_shift,
              flip,
              bckg_class,
              std::string("b_bird!"),
              n_threads) ;
em.classify() ;
Matrix4D<double> prob  = em.get_post_prob() ;
Matrix3D<double> model = em.get_sequence_models() ;


// procedure to extend the models
// 1 modify from/to range
int   n_ext       = 1000  ; // add 500 columns of 1bp
int   n_ext_right = n_ext/2 ;
int   n_ext_left  = n_ext - n_ext_right ;
from -= ext_left ;
to   += ext_right ;

// 2 create the extended matrix
CorrelationMatrixCreator mc2(bed_path,
                             bam_path,
                             bai_path,
                             from,
                             to) ;
Matrix2D<int> data_ext = mc2.create_matrix() ;

// 3 compute the models
SequenceModelComputer smc(data_ext,
                          prob,
                          bckg_class,
                          n_threads) ;
Matrix2D<double> model_ext = smc.get_model() ;

// 4 display the extended models
std::cout << model_ext << std::endl ;

```

The same could be performed using the autonomous application SequenceModelExtender that implements this procedure.

```
SequenceModelExtender --bed /path/to/file.bed --fasta /path/to/file.fasta --prob /path/to/prob.mat4d --from -1000 --to 1000 --binSize 1 --ext 1000 --thread 4
```

For help, enter :

```
SequenceModelExtender -h
```

## 4.3 Extracting data assigned to a class

![](images/extraction.png?raw=true)
**Class data extraction procedure**
Toy example showing the data extraction procedure by retaking the previous example.
On the left, the procedure for read density data is shown. A data matrix is created using the BED, BAM and BAI files. The from/to range determine the final number of bins. Because shifting trims the number of columns, the initial matrix is extended to account for this. Then, instead of creating a weighted aggregation as to compute class models, a per region (row) weighted aggregation is performed. For a given region, all its slices are aggregated using the probabilities as weights. This results in the creation of one vector per region. Eventually, each of these vectors corresponds to one row (region) in the final 2D matrix.
On the right, the procedure for DNA sequences is shown. It works exactly the same as for read density data with the exception that it creates a 3D matrix containing one consensus sequences per row (region). Indeed, the per region aggregation procedure creates a LPM instead of a vector.
On the edges, the C++ classes that implemente these procedures, their inputs and their outputs are displayed for read density and DNA sequence data.


The aim of the manipulation described in this section conceptually corresponds to creating a matrix containing only the rows (region) assigned to a given class X. This is of interest to run downstream analysis on this specific data subset (further clustering for instance). For hard clustering algorithms, this can be done quite easily. It only requires to select the regions (rows) that have been assigned to class X. However, for soft partitioning, things are different and each and every region is assigned to all classes, with varying probabilities. This procedure is displayed in the figure above.

### 4.3.1 Extracting a read density class matrix

The following snippet of code indicates how to extract the data assigned to each class of a given partition.

```
#include <iostream>
#include <string>
#include "GenomicTools/SequenceMatrixCreator.hpp"
#include "GenomicTools/ClassReadDataCreator.hpp"
#include "Clustering/EMRead.hpp"
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix4D.hpp"

// create data matrix
std::string file_bed("/some/bed/file.bed") ;
std::string file_bam("/some/bam/file.bam") ;
std::string file_bai("/some/bai/file.bai") ;

int from     = -1000 ;
int to       =  1000 ;
int bin_size =    10 ;

std::string method = "read" ;

CorrelationMatrixCreator mc(file_bed,
                            file_bam,
                            file_bai,
                            from,
                             to,
                             bin_size,
                            method) ;
Matrix2D<int> data = mc.create_matrix() ;

// partition the data
size_t n_class    = 4 ;
size_t n_iter     = 20 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;
std::string file_prob = "prob.mat4d" ;

EMRead em(data,
          n_class,
          n_iter,
          n_shift,
          flip,
          bckg_class,
          "looking for a seed...",
          n_threads) ;
em.classify() ;
Matrix4D<double> prob = em.get_post_prob() ;
prob.save(file_prob) ;

// extract each class and displays it
for(size_t class=1; class<=n_class; class++)
{	ClassReadDataCreator crc(file_bed,
 	                         file_bam,
 	                         file_bai,
 	                         file_prob,
 	                         from,
 	                         to,
 	                         bin_size,
 	                         class,
 	                         method) ;
	Matrix2D<double> class_mat = crc.create_matrix()
	std::cout << class_mat << std::endl ;
}

```
The same can be done using the autonomous application ClassReadDataCreator.

```
ClassReadDataCreator --bed /some/bed/file.bed --bam /some/bed/file.bam --bai /some/bed/file.bai --prob /some/bed/prob.mat4d --from -1000 --to 1000 --binSize 10 --k <class> --method read
```

For help, enter :
```
ClassReadDataCreator -h
```


### 4.3.2 Extracting a DNA sequence class matrix

The following snippet of code indicates how to extract the data assigned to each class of a given partition.

```
#include <iostream>
#include <string>
#include "GenomicTools/SequenceMatrixCreator.hpp"
#include "GenomicTools/ClassSequenceDataCreator.hpp"
#include "Clustering/EMSequence.hpp"
#include "Matrix/Matrix2D.hpp"
#include "Matrix/Matrix3D.hpp"
#include "Matrix/Matrix4D.hpp"

// create data matrix
std::string file_bed("/some/bed/file.bed") ;
std::string file_seq("/some/bam/file.fasta") ;

int from     = -1000 ;
int to       =  1000 ;

SequenceMatrixCreator mc(file_bed,
                         file_seq,
                         from,
                         to) ;
Matrix2D<int> data = mc.create_matrix() ;

// partition the data
size_t n_class    = 4 ;
size_t n_iter     = 20 ;
size_t n_shift    = 21 ;
bool   flip       = true ;
bool   bckg_class = false ;
size_t n_threads  = 4 ;
std::string file_prob = "prob.mat4d" ;

EMSequence em(data,
              n_class,
              n_iter,
              n_shift,
              flip,
              bckg_class,
              "your genome is full of mistakes",
              n_threads) ;
em.classify() ;
Matrix4D<double> prob = em.get_post_prob() ;
prob.save(file_prob) ;

// extract each class and displays it
for(size_t class=1; class<=n_class; class++)
{	ClassSequenceDataCreator csc(file_bed,
 	                             file_seq,
 	                             file_prob,
 	                             from,
 	                             to,
 	                             class) ;
	Matrix3D<double> class_mat = csc.create_matrix()
	std::cout << class_mat << std::endl ;
}

```
The same can be done using the autonomous application ClassSequenceDataCreator.

```
ClassSequenceDataCreator --bed /some/bed/file.bed --fasta /some/bed/file.fasta --prob /some/bed/prob.mat4d --from -1000 --to 1000 --k <class>
```

For help, enter :
```
ClassSequenceDataCreator -h
```


# 5 Dependencies

The code present on this repository relies on the following external libraries :

1. boost v1.4.1 or higher (https://www.boost.org/)
2. UnitTest++ v2 (https://github.com/unittest-cpp/unittest-cpp)
3. SeqAn v2.4.0 or higher (https://seqan.readthedocs.io/en/master/)

for command line parsing, for some string and mathematial operations (boost), for unit testing (UnitTest++) and for I/O operations on BAM, BED and FASTA files (SeqAn). These libraries are not provided with this project and need to be installed. More information about their installation is provided in the following section.


# 6 Install and compile

For linux systems, scripts installing the above mentionned libraries are provided. In the scripts/ directory, 3 shell scripts allow to install these libraries. It is highly recommanded to use them. By doing so, the libraries will be installed in a sub-directory lib/ at the root of the project directory. To do so, open a shell script and type :

```
cd <path_to_project_root>

# install boost
scripts/install_libboost.sh

#install seqan
scripts/install_libSeqAn.sh

#install unittest++
scripts/install_libUnitTest++.sh
```

Then, run 

```
cd <path_to_project_root>
scripts/config.sh

```

It will replace all occurences of "<path_to_project_root>" by the path to the project root, for instance "/local/user/EMtools" in ./build.sh, ./CMakeLists.txt and src/CMakeLists.txt files. Once this is done, the project can be compiled using :

```
cd <path_to_project_root>

# compile
./build.sh
```

At the end, the C++ static libaries of each component of the project will be found in the lib/ directory. The application executables will be stored in the bin/ directory.


# 6 Authors

* **Romain Groux**


# 7 License

This project is licensed under the GNU General Public License v3 - see the [LICENSE.md](LICENSE.md) file for details


# 8 Acknowledgments

* Vincent Gardeux
* Philipp Bucher
* Giovanna Ambosini


