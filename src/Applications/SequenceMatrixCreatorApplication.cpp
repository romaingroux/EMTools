
#include <SequenceMatrixCreatorApplication.hpp>
#include <SequenceMatrixCreator.hpp>

#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <stdexcept>                   // std::invalid_argument


namespace po = boost::program_options ;

// the valid values for --method option
std::string method_read            = "read" ;
std::string method_read_atac       = "read_atac" ;
std::string method_fragment        = "fragment" ;
std::string method_fragment_center = "fragment_center" ;


CorrelationMatrixCreatorApplication::CorrelationMatrixCreatorApplication(int argn, char** argv)
    : file_bed(""), file_fasta(""), from(0), to(0), runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int CorrelationMatrixCreatorApplication::run()
{   if(this->runnable)
    {   SequenceMatrixCreator mc(this->file_bed,
                                 this->file_fasta,
                                 this->from,
                                 this->to) ;

        std::cout << mc.create_matrix() << std::endl ;
        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE ; }
}

void CorrelationMatrixCreatorApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "SequenceMatrixCreator is an application that creates a\n"
                                   "sequence matrix (1 nucleotide per matrix cell in int format)\n"
                                   "from a BED file and a fasta file and returnes it through stdout.\n"
                                   "The sequence centers are defined as the center position of each\n"
                                   "region contained in the bed file. The corresponding single bp\n"
                                   "regions are then extended using the from/to parameters on each side.\n"
                                   "The corresponding sequences are then extracted from the fasta file.\n"
                                   "The sequences are found by searching the sequences in the fasta file\n. "
                                   "Thus, the sequence IDs in the BED file are expected to match EXACTLY\n"
                                   "the sequence headers in the fasta file.\n"
                                   "The DNA character to integer conversion code is the following :\n"
                                   "A->0, C->1, G->2, T->3, N->4.\n"
                                   "The matrix is a 2D matrix which dimensions are :\n"
                                   "1) number of regions\n"
                                   "2) length of region (<to> - <from> + 1)\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_bed_msg      = "The path to the BED file containing the references.";
    std::string opt_fasta_msg    = "The path to the fasta file containing the target sequences. The file "
                                   "extension must be .fa or .fasta.";
    std::string opt_from_msg     = "The upstream limit - in relative coordinate - of the region to build "
                                   "around each reference center." ;
    std::string opt_to_msg       = "The downstream limit - in relative coordinate - of the region to build "
                                   "around each reference center." ;

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    std::string method(method_read) ;

    desc.add_options()
                ("help,h", opt_help_msg.c_str())

                ("bed",    po::value<std::string>(&(this->file_bed)), opt_bed_msg.c_str())
                ("fasta",  po::value<std::string>(&(this->file_fasta)), opt_fasta_msg.c_str())

                ("from",   po::value<int>(&(this->from)),             opt_from_msg.c_str())
                ("to",     po::value<int>(&(this->to)),               opt_to_msg.c_str()) ;

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
    if(this->file_bed == "" and (not help))
    {   std::string msg("Error! No BED file was given (--bed)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->file_fasta == "" and (not help))
    {   std::string msg("Error! No fasta file was given (--fasta)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->from == 0 and this->to == 0 and (not help))
    {   std::string msg("Error! No range given (--from and --to)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->from >= this->to and (not help))
    {   std::string msg("Error! from shoud be smaller than to (--from and --to)!") ;
        throw std::invalid_argument(msg) ;
    }

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
{   CorrelationMatrixCreatorApplication app(argn, argv) ;
    return app.run() ;
}

