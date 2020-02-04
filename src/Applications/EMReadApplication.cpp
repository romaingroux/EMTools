
#include <EMReadApplication.hpp>
#include <EMRead.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <utility>                     // std::move()
#include <stdexcept>                   // std::invalid_argument
#include <boost/program_options.hpp>

#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <matrix_utility.hpp>   // filter()


namespace po = boost::program_options ;


EMReadApplication::EMReadApplication(int argn, char** argv)
    : file_read(""), file_filter(""), file_out(""),
      n_class(0), n_iter(0), n_shift(0), flip(false), bckg_class(false),
      n_threads(0), seed(""), runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int EMReadApplication::run()
{   if(this->runnable)
    {
        // read data
        Matrix2D<int> data(this->file_read) ;

        // filter out some rows if needed
        std::vector<size_t> filter ;
        if(this->file_filter != "")
        {   // it is a column vector, easier to use the Matrix2D interface
            // to read it rather than coding a function for :)
            filter = Matrix2D<size_t>(this->file_filter).get_data() ;
            std::sort(filter.begin(), filter.end()) ;
            data = filter_rows(filter, data) ;
        }

        EMRead em(std::move(data),
                  this->n_class,
                  this->n_iter,
                  this->n_shift,
                  this->flip,
                  this->bckg_class,
                  this->seed,
                  this->n_threads) ;
        em.classify() ;
        em.get_post_prob().save(this->file_out) ;
        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE ; }
}

void EMReadApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "EMRead is a probabilistic partitioning algorithm that \n"
                                   "sofetly assigns genomic regions to classes given the shape \n"
                                   "of the read/fragment count density over the region. The\n"
                                   "assignment probabilities are written in binary format as a"
                                   "4D matrix of dimensions. Its dimensions are :\n"
                                   "1) number of regions\n"
                                   "2) number of classes \n"
                                   "3) number of shift states\n"
                                   "4) number of flip states\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_thread_msg   = "The number of threads dedicated to parallelize the computations,\n "
                                   "by default 0 (no parallelization)." ;
    std::string opt_read_msg     = "The path to the file containing the read density data" ;
    std::string opt_filter_msg   = "Optional. The path to a single column text file containing the 0-based\n"
                                   "indices of rows to filter out in the data." ;
    std::string opt_file_out_msg = "A path to a file in which the assignment probabilities will be saved\n"
                                   "in binary format." ;
    std::string opt_iter_msg     = "The number of iterations." ;
    std::string opt_class_msg    = "The number of classes to find." ;
    std::string opt_shift_msg    = "Enables this number of column of shifting "
                                   "freedom to realign the data. By default, shifting is "
                                   "disabled (equivalent to --shift 1)." ;
    std::string opt_flip_msg     = "Enables flipping to realign the data.";
    std::string opt_bckg_msg     = "Adds a class to model the signal background. This class\n"
                                   "contains the mean number of read at each position and\n"
                                   "is never updated." ;
    std::string opt_seed_msg     = "A value to seed the random number generator.";

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    std::string seeding_tmp ;

    desc.add_options()
                ("help,h",   opt_help_msg.c_str())

                ("read",     po::value<std::string>(&(this->file_read)),   opt_read_msg.c_str())
                ("filter",   po::value<std::string>(&(this->file_filter)), opt_filter_msg.c_str())
                ("out",      po::value<std::string>(&(this->file_out)),    opt_file_out_msg.c_str())

                ("iter,i",   po::value<size_t>(&(this->n_iter)),           opt_iter_msg.c_str())
                ("class,c",  po::value<size_t>(&(this->n_class)),          opt_class_msg.c_str())
                ("shift,s",  po::value<size_t>(&(this->n_shift)),          opt_shift_msg.c_str())
                ("flip",     opt_flip_msg.c_str())
                ("bgclass",  opt_bckg_msg.c_str())

                ("seed",     po::value<std::string>(&(this->seed)),        opt_seed_msg.c_str())
                ("thread",   po::value<std::size_t>(&(this->n_threads)),   opt_thread_msg.c_str()) ;

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
    if(this->file_read == "" and
       (not help))
    {   std::string msg("Error! No data were given (--read)!") ;
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
{   EMReadApplication app(argn, argv) ;
    return app.run() ;
}

