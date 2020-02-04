
#include <ClassSequenceDataCreatorApplication.hpp>
#include <ClassSequenceDataCreator.hpp>

#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>                   // std::invalid_argument

#include <Matrix4D.hpp>

namespace po = boost::program_options ;

ClassSequenceDataCreatorApplication::ClassSequenceDataCreatorApplication(int argn, char** argv)
    : file_bed(""), file_fasta(""), file_prob(""), file_out(""),
      from(0), to(0), class_k(0),
      runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int ClassSequenceDataCreatorApplication::run()
{   if(this->runnable)
    {
        ClassSequenceDataCreator csc(this->file_bed,
                                     this->file_fasta,
                                     this->file_prob,
                                     this->from,
                                     this->to,
                                     this->class_k) ;

        // display integer matrix for given class
        csc.create_matrix().save(this->file_out) ;
        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE ; }
}

void ClassSequenceDataCreatorApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "ClassSequenceDataCreator is an autonomous application that extracts the data\n"
                                   "that have been assigned to a given class K.\n"
                                   "\n"
                                   "Given posterior probabilities and a sequence matrix, the corresponding\n"
                                   "class models can be computed. They are the weighted aggregations of the\n"
                                   "DNA sequences assigned to each given class. However, because DNA sequences\n"
                                   "cannot be summed, the aggregation are represented as probability matrices\n"
                                   "or consensus sequence (A+C is represented as 50%A, 50%C, 0%G, 0%T). Instead\n"
                                   "of this, this program creates the unfolded matrix that, if summed over the\n"
                                   "columns, gives the model of class K.\n"
                                   "\n"
                                   "For a hard clustering methods, this procedure would simply correspond to the\n"
                                   "creation of a matrix of dimensions N'xL where N'<=N is the number sequences\n"
                                   "assigned to class K among the N overall sequences and L the length of\n"
                                   "the each sequence.\n"
                                   "\n"
                                   "In the case of a soft clustering methods, this procedure creates a 3D matrix of\n"
                                   "dimensions NxL'x4. This matrix contains N probability matrices, each one of \n"
                                   "dimensions L'x4 where L'=L-S+1, 4 corresponds to A, C, G, T and S is the\n"
                                   "shifting freedom allowed during the classification. The resulting matrix\n"
                                   "contains as many rows as the starting matrix because in soft clustering, all\n"
                                   "sequences (rows) are assigned to all classes\n"
                                   "\n"
                                   "To construct a final matrix M3 of dimensions NxL3 where L3 covers a given\n"
                                   "range <from>/<to>, the original matrix M1 of dimensions NxL is computed and\n"
                                   "extended into a matrix M2 NxL2 with L2>=L1. The final M3 of dimensions NxL\n"
                                   "is eventually computed, for class K, using the given posterior probabilities.\n"
                                   "A row of the final matrix M3 is the weighted average of each of the S\n"
                                   "possibles slices of the corresponding row in M2, represented as a probability\n"
                                   "matrix. The weights used are the probabilities with which this row was assigned\n"
                                   "to class K, for each of the S shift states, in each flip state.\n"
                                   "\n"
                                   "The original matrix M1 that was partitionned with shifting freedom S is\n"
                                   "generated using the BED and fasta files that were originally used to\n"
                                   "create it.\n"
                                   "The posterior probabilities should be a 4D matrix in binary format, with\n"
                                   "dimensions :\n"
                                   "1) number of sequences\n"
                                   "2) number of classes\n"
                                   "3) number of shift states\n"
                                   "4) number of flip states\n"
                                   "The results is returned as a 3D binary matrix of dimensions :\n"
                                   "1) number of sequences\n"
                                   "2) length of the sequences, as defined by the <from>/<to> range\n"
                                   "3= 4 for A, C, G, T\n\n" ;

    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_bed_msg      = "The path to the BED file containing the references.";
    std::string opt_fasta_msg    = "The path to the fasta file containing the sequences.";
    std::string opt_prob_msg     = "The path to the file containing the assignment probabilities\n"
                                   "(the partition)." ;
    std::string opt_out_msg      = "The path to the file in which the results will be written (3D matrix\n"
                                   "in binary format)." ;
    std::string opt_from_msg     = "The most upstream position - in relative coordinate - of the regions\n"
                                   "in the original matrix." ;
    std::string opt_to_msg       = "The most downstream position - in relative coordinate - of the regions\n"
                                   "in the original matrix." ;
    std::string opt_classk_msg   = "The index (1-based) of the class of interest." ;

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    desc.add_options()
                ("help,h",  opt_help_msg.c_str())

                ("bed",     po::value<std::string>(&(this->file_bed)),   opt_bed_msg.c_str())
                ("fasta",   po::value<std::string>(&(this->file_fasta)), opt_fasta_msg.c_str())
                ("prob",    po::value<std::string>(&(this->file_prob)),  opt_prob_msg.c_str())

                ("out",     po::value<std::string>(&(this->file_out)),   opt_out_msg.c_str())

                ("from",    po::value<int>(&(this->from)),               opt_from_msg.c_str())
                ("to",      po::value<int>(&(this->to)),                 opt_to_msg.c_str())
                ("k",       po::value<size_t>(&(this->class_k)),         opt_classk_msg.c_str()) ;

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
    else if(this->file_prob == "" and (not help))
    {   std::string msg("Error! No probability (partition) file was given (--prob)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->file_out == "" and (not help))
    {   std::string msg("Error! No output file was given (--out)!") ;
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
    else if(this->class_k == 0 and (not help))
    {   std::string msg("Error! no class of interest was given (--k)!") ;
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
{   ClassSequenceDataCreatorApplication app(argn, argv) ;
    return app.run() ;
}

