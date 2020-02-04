
#include <EMConsensusSequenceApplication.hpp>
#include <EMConsensusSequence.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <utility>                     // std::move()
#include <stdexcept>                   // std::invalid_argument
#include <boost/program_options.hpp>

#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <matrix_utility.hpp>
#include <KmersStatistics.hpp>  // kmer::compute_kmer_pvalue()
#include <sorting_utility.hpp>  // order()


namespace po = boost::program_options ;


EMConsensusSequenceApplication::EMConsensusSequenceApplication(int argn, char** argv)
    : file_consseq(""), file_filter(""), file_out(""),
      n_class(0), n_iter(0), n_shift(0), flip(false), bckg_class(false),
      n_threads(0), seed(""), runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int EMConsensusSequenceApplication::run()
{   if(this->runnable)
    {   EMConsensusSequence* em(nullptr) ;

        // row filter
        std::vector<size_t> filter ;
        if(this->file_filter != "")
        {   // it is a column vector, easier to use the Matrix2D interface
            // to read it rather than coding a function for :)
            filter = Matrix2D<size_t>(this->file_filter).get_data() ;
            std::sort(filter.begin(), filter.end()) ;
        }

        // data
        Matrix3D<double> data ;
        data.load(this->file_consseq) ;
        // filter out some rows if needed
        if(filter.size())
        {   data = filter_rows(filter, data) ; }

        // seeds motifs randomly
        if(this->seed != "")
        {   em = new EMConsensusSequence(std::move(data),
                                         this->n_class,
                                         this->n_iter,
                                         this->n_shift,
                                         this->flip,
                                         this->bckg_class,
                                         this->seed,
                                         this->n_threads) ;
        }
        // seeds from enriched kmers
        else
        {   size_t model_ncol = data.get_dim()[1] - this->n_shift + 1 ;

            Matrix3D<double> model = this->init_model_kmer(model_ncol,
                                                           data) ;
            em = new EMConsensusSequence(std::move(data),
                                         std::move(model),
                                         this->n_iter,
                                         this->flip,
                                         this->bckg_class,
                                         this->n_threads) ;
        }

        // classify
        em->classify() ;
        em->get_post_prob().save(this->file_out) ;

        // clean
        delete em ;
        em = nullptr ;

        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE ; }
}

void EMConsensusSequenceApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "EMConsensusSequence is a probabilistic partitioning algorithm that \n"
                                   "sofetly assigns DNA consensus sequences to classes given their motif\n"
                                   "content.\n"
                                   "A DNA consensus sequence is a probability matrix that contains the\n"
                                   "probability of each base at each position. The input is thus a 3D\n"
                                   "matrix that contains several consensus matrices (that are 2D\n"
                                   "matrices\n"
                                   "In essence, this application does exactly the same as EMSequence,\n"
                                   "with the exception that it takes consensus sequences instead of\n"
                                   "regular DNA sequences. The consensus sequence matrix should be\n"
                                   "stored in binary format. Its dimensions should be :\n"
                                   "1) the number of sequences\n"
                                   "2) the length of the sequences\n"
                                   "3) 4 for A, C, G, T\n"
                                   "The probabilities for a given posision (2nd dimension) in a consensus\n"
                                   "sequence (1st dimensions) must sum to 1."
                                   "The assignment probabilities are written in binary format as a 4D "
                                   "matrix which dimensions are :\n"
                                   "1) the number of sequences\n"
                                   "2) the number of classes\n"
                                   "2) the number of shift states\n"
                                   "3) the number of flip states\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_thread_msg   = "The number of threads dedicated to parallelize the computations,\n "
                                   "by default 0 (no parallelization)." ;
    std::string opt_consseq_msg  = "The path to the file containing the consensus sequences" ;
    std::string opt_filter_msg   = "Optional. The path to a single column text file containing the 0-based\n"
                                   "indices of rows to filter out in the data." ;
    std::string opt_file_out_msg = "A path to a file in which the assignment probabilities will be saved\n"
                                   "in binary format." ;
    std::string opt_iter_msg     = "The number of iterations." ;
    std::string opt_class_msg    = "The number of classes to find." ;
    std::string opt_shift_msg    = "Enables this number of column of shifting freedom to realign\n"
                                   "the data. By default, shifting is disabled (equivalent to\n"
                                   "--shift 1)." ;
    std::string opt_flip_msg     = "Enables flipping to realign the data.";
    std::string opt_bckg_msg     = "Adds a class to model the sequence background. This class\n"
                                   "contains the sequence background probabilities at each position\n"
                                   "and is never updated." ;
    std::string opt_seed_msg     = "A value to seed the random number generator.";

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    std::string seeding_tmp ;

    desc.add_options()
                ("help,h",  opt_help_msg.c_str())

                ("consseq", po::value<std::string>(&(this->file_consseq)), opt_consseq_msg.c_str())
                ("filter",   po::value<std::string>(&(this->file_filter)), opt_filter_msg.c_str())

                ("out",     po::value<std::string>(&(this->file_out)),     opt_file_out_msg.c_str())

                ("iter,i",  po::value<size_t>(&(this->n_iter)),            opt_iter_msg.c_str())
                ("class,c", po::value<size_t>(&(this->n_class)),           opt_class_msg.c_str())
                ("shift,s", po::value<size_t>(&(this->n_shift)),           opt_shift_msg.c_str())
                ("flip",    opt_flip_msg.c_str())
                ("bgclass", opt_bckg_msg.c_str())

                ("seed",    po::value<std::string>(&(this->seed)),         opt_seed_msg.c_str())
                ("thread",  po::value<std::size_t>(&(this->n_threads)),    opt_thread_msg.c_str()) ;

    // parse
    try
    {   po::store(po::parse_command_line(argn, argv, desc), vm) ;
        po::notify(vm) ;
    }
    catch(std::invalid_argument& e)
    {   std::string msg = std::string("Error! Invalid option given!\n") + std::string(e.what()) ;
        throw std::invalid_argument(msg) ;
    }
    catch(...)
    {	throw std::invalid_argument("An unknown error occured while parsing the options") ; }

    bool help = vm.count("help") ;

    // checks unproper option settings
    if(this->file_consseq == "" and
       (not help))
    {   std::string msg("Error! No data were given (--seq)!") ;
        throw std::invalid_argument(msg) ;
    }
    if(this->file_out == "" and
       (not help))
    {   std::string msg("Error! No output file given (--out)!") ;
        throw std::invalid_argument(msg) ;
    }

    // no iter given -> 1 iter
    if(this->n_iter == 0)
    {   this->n_iter = 1 ; }
    // no shift class given -> 1 class
    if(this->n_class == 0)
    {   this->n_class = 1 ; }
    // no shift given, value of 1 -> no shift
    if(this->n_shift == 0)
    {   this->n_shift = 1 ; }
    // set flip
    if(vm.count("flip"))
    {   this->flip  = true ; }
    // set background class
    if(vm.count("bgclass"))
    {   this->bckg_class  = true ; }
    // help invoked, run() cannot be invoked
    if(help)
    {   std::cout << desc << std::endl ;
        this->runnable = false ;
        return ;
    }
    // everything fine, run() can be called
    else
    {   this->runnable = true ;
        return ;
    }
}

Matrix3D<double> EMConsensusSequenceApplication::init_model_kmer(size_t l_model,
                                                                 const Matrix3D<double>& data) const
{
    // leave space for 2N's on each side
    size_t l_kmer  = l_model ;
    size_t n_n     = 0 ; // so far, 0 N's added
    if(l_model > 4)
    {   n_n     = 2 ; // 2 N's on each side
        l_kmer -= (2*n_n) ;
    }

    // compute the pvalue associated to each kmer
    auto kmers_pvalues = kmers::compute_kmer_pvalue(data, l_kmer) ;

    // sort kmers by ascending pvalue
    std::vector<size_t> index = order(kmers_pvalues.second, true) ;
    // get most significant
    std::vector<std::string> kmers(this->n_class) ;
    for(size_t i=0; i<this->n_class; i++)
    {   size_t idx = index[i] ;
        kmers[i] = kmers_pvalues.first[idx] ;
        std::cerr << kmers_pvalues.first[idx] << "  " << kmers_pvalues.second[idx] << std::endl ;
    }
    // turn to motifs
    double p_base  = 0.7 ;  // the prob of the base matching these of the kmer
    double p_nbase = 0.1 ;  // the prob of the bases not matching these of the kmer
    double p_n     = 0.25 ; // the prob of N
    // only N's for now
    Matrix3D<double> model(this->n_class,
                                l_model,
                                4,
                                p_n) ;
    for(size_t i=0; i<kmers.size(); i++)
    {   for(size_t j_kmer=0, j_model=n_n; j_kmer<l_kmer; j_kmer++)
        {   // A
            if(kmers[i][j_kmer] == '0')
            {   model(i,j_model, 0) = p_base ;
                model(i,j_model, 1) = p_nbase ;
                model(i,j_model, 2) = p_nbase ;
                model(i,j_model, 3) = p_nbase ;
            }
            // C
            else if(kmers[i][j_kmer] == '1')
            {   model(i,j_model, 0) = p_nbase ;
                model(i,j_model, 1) = p_base ;
                model(i,j_model, 2) = p_nbase ;
                model(i,j_model, 3) = p_nbase ;
            }
            // G
            else if(kmers[i][j_kmer] == '2')
            {   model(i,j_model, 0) = p_nbase ;
                model(i,j_model, 1) = p_nbase ;
                model(i,j_model, 2) = p_base ;
                model(i,j_model, 3) = p_nbase ;
            }
            // T
            else if(kmers[i][j_kmer] == '3')
            {   model(i,j_model, 0) = p_nbase ;
                model(i,j_model, 1) = p_nbase ;
                model(i,j_model, 2) = p_nbase ;
                model(i,j_model, 3) = p_base ;
            }
        }
    }
    return model ;
}

int main(int argn, char** argv)
{   EMConsensusSequenceApplication app(argn, argv) ;
    return app.run() ;
}

