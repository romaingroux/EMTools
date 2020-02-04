
#include <EMJointApplication.hpp>
#include <EMJoint.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <utility>                     // std::move()
#include <stdexcept>                   // std::invalid_argument
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>  // boost::split()

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <matrix_utility.hpp>   // filter()

namespace po = boost::program_options ;


EMJointApplication::EMJointApplication(int argn, char** argv)
    : files_read(""), file_sequence(""), file_conssequence(""),
      file_filter(""), file_out(""),
      n_class(0), n_iter(0), n_shift(0), flip(false), bckg_class(false),
      n_threads(0), seed(""), runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int EMJointApplication::run()
{   if(this->runnable)
    {   // read data
        std::vector<std::string> read_paths ;
        boost::split(read_paths, this->files_read, [](char c){return c == ',';}) ;

        // row filter
        std::vector<size_t> filter ;
        if(this->file_filter != "")
        {   // it is a column vector, easier to use the Matrix2D interface
            // to read it rather than coding a function for :)
            filter = Matrix2D<size_t>(this->file_filter).get_data() ;
            std::sort(filter.begin(), filter.end()) ;
        }

        std::vector<Matrix2D<int>> data_read ;
        for(const auto& path : read_paths)
        {   if(path == "")
            {   continue ; }
            Matrix2D<int> data(path) ;
            // filter out some rows if needed
            if(filter.size())
            {   data = filter_rows(filter, data) ; }
            data_read.push_back(std::move(data)) ;
        }
        // read data only
        EMJoint* em = nullptr ;
        if(this->file_sequence     == "" and
           this->file_conssequence == "")
        {   em = new EMJoint(std::move(data_read),
                             this->n_class,
                             this->n_iter,
                             this->n_shift,
                             this->flip,
                             this->bckg_class,
                             this->seed,
                             this->n_threads) ;
        }
        // read and sequence data
        else if(this->file_sequence     != "" and
                this->file_conssequence == "")
        {   Matrix2D<int> data_seq(this->file_sequence) ;
            // filter out some rows if needed
            if(filter.size())
            {   data_seq = filter_rows(filter, data_seq) ; }

            em = new EMJoint(std::move(data_read),
                             std::move(data_seq),
                             this->n_class,
                             this->n_iter,
                             this->n_shift,
                             this->flip,
                             this->bckg_class,
                             this->seed,
                             this->n_threads) ;
        }
        // read and consensus sequence data
        else if(this->file_sequence     == "" and
                this->file_conssequence != "")
        {   Matrix3D<double> data_consseq ;
            data_consseq.load(this->file_conssequence) ;
            // filter out some rows if needed
            if(filter.size())
            {   data_consseq = filter_rows(filter, data_consseq) ; }
            em = new EMJoint(std::move(data_read),
                             std::move(data_consseq),
                             this->n_class,
                             this->n_iter,
                             this->n_shift,
                             this->flip,
                             this->bckg_class,
                             this->seed,
                             this->n_threads) ;
        }
        em->classify() ;
        em->get_post_prob().save(this->file_out) ;
        delete em ;
        em = nullptr ;
        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE ; }
}

void EMJointApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "EMJoint is a probabilistic partitioning algorithm that \n"
                                   "sofetly assigns genomic regions to classes given 1) the shapes \n"
                                   "of the read/fragment count densities over the regions and\n"
                                   "2) the region DNA sequences.\n "
                                   "The assignment probabilities are written in a binary format as\n"
                                   "a 4D matrix of dimensions number of regions x number of\n"
                                   "classes x number of shift states x number of flip states\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_thread_msg   = "The number of threads dedicated to parallelize the computations, \n"
                                   "by default 0 (no parallelization)." ;
    std::string opt_read_msg     = "A coma separated list of paths to the file containing the \n"
                                   "read density data. At least one path is needed." ;
    std::string opt_seq_msg      = "The path to the file containing the sequence data. If no path is \n"
                                   "given, the classification is only cares about the read density shapes." ;
    std::string opt_consseq_msg  = "The path to the file containing the consensus sequence data. If no path is \n"
                                   "given, the classification only cares about the read density shapes." ;
    std::string opt_filter_msg   = "Optional. The path to a single column text file containing the 0-based\n"
                                   "indices of rows to filter out in the data." ;
    std::string opt_file_out_msg = "A path to a file in which the assignment probabilities will be saved\n"
                                   "in binary format." ;
    std::string opt_iter_msg     = "The number of iterations." ;
    std::string opt_class_msg    = "The number of classes to find." ;
    std::string opt_shift_msg    = "Enables this number of column of shifting "
                                   "freedom. By default, shifting is "
                                   "disabled (equivalent to --shift 1)." ;
    std::string opt_flip_msg     = "Enables flipping.";
    std::string opt_bckg_msg     = "Adds a class to model the sequence and the read signal background. This\n"
                                   "class contains sequence background probabilies (for the sequence model)\n"
                                   "and the mean number of reads (for the read count models) at each position\n"
                                   "and is never updated." ;
    std::string opt_seed_msg     = "A value to seed the random number generator.";

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    desc.add_options()
                ("help,h",   opt_help_msg.c_str())

                ("read",     po::value<std::string>(&(this->files_read)),        opt_read_msg.c_str())
                ("seq",      po::value<std::string>(&(this->file_sequence)),     opt_seq_msg.c_str())
                ("consseq",  po::value<std::string>(&(this->file_conssequence)), opt_seq_msg.c_str())
                ("filter",   po::value<std::string>(&(this->file_filter)),       opt_filter_msg.c_str())
                ("out",      po::value<std::string>(&(this->file_out)),          opt_file_out_msg.c_str())

                ("iter,i",   po::value<size_t>(&(this->n_iter)),                 opt_iter_msg.c_str())
                ("class,c",  po::value<size_t>(&(this->n_class)),                opt_class_msg.c_str())
                ("shift,s",  po::value<size_t>(&(this->n_shift)),                opt_shift_msg.c_str())
                ("flip",     opt_flip_msg.c_str())
                ("bgclass",  opt_bckg_msg.c_str())

                ("seed",     po::value<std::string>(&(this->seed)),              opt_seed_msg.c_str())
                ("thread",   po::value<std::size_t>(&(this->n_threads)),         opt_thread_msg.c_str()) ;

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
    if(this->files_read == "" and
       this->file_sequence == "" and
       this->file_conssequence == "" and
       (not help))
    {   std::string msg("Error! No data were given (--read --seq --consseq)!") ;
        throw std::invalid_argument(msg) ;
    }
    if(this->files_read == "" and
       (not help))
    {   std::string msg("Error! No read density data were given (--read)!") ;
        throw std::invalid_argument(msg) ;
    }
    if(this->file_sequence != "" and
       this->file_conssequence != "")
    {   std::string msg("Error! --seq and --consseq are mutually exclusive!") ;
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


int main(int argn, char** argv)
{   EMJointApplication app(argn, argv) ;
    return app.run() ;
}

