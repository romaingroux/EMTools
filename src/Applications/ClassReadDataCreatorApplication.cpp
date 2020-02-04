
#include <ClassReadDataCreatorApplication.hpp>
#include <ClassReadDataCreator.hpp>

#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>                   // std::invalid_argument

#include <Matrix4D.hpp>

namespace po = boost::program_options ;

// the valid values for --method option
std::string method_read            = "read" ;
std::string method_read_atac       = "read_atac" ;
std::string method_fragment        = "fragment" ;
std::string method_fragment_center = "fragment_center" ;


ClassReadDataCreatorApplication::ClassReadDataCreatorApplication(int argn, char** argv)
    : file_bed(""), file_bam(""), file_bai(""), file_prob(""),
      from(0), to(0), bin_size(0), class_k(0),
      method(CorrelationMatrixCreator::FRAGMENT), runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int ClassReadDataCreatorApplication::run()
{   if(this->runnable)
    {   ClassReadDataCreator crc(this->file_bed,
                                 this->file_bam,
                                 this->file_bai,
                                 this->file_prob,
                                 this->from,
                                 this->to,
                                 this->bin_size,
                                 this->class_k,
                                 this->method) ;

        // display integer matrix for given class
        std::cout << crc.create_matrix() << std::endl ;
        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE ; }
}

void ClassReadDataCreatorApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "ClassReadDataCreator is an autonomous application that extracts the data\n"
                                   "that have been assigned to a given class K.\n"
                                   "\n"
                                   "Given posterior probabilities and a read density matrix, the corresponding\n"
                                   "class models can be computed. They are the weighted aggregations of the\n"
                                   "data assigned to each given class. Instead of this, this program creates\n"
                                   "the unfolded matrix that, if summed over the columns, gives the model of\n"
                                   "class K.\n"
                                   "\n"
                                   "For a hard clustering methods, this procedure would simply correspond to the\n"
                                   "creation of a matrix of dimensions N'xL where N'<=N is the number regions\n"
                                   "assigned to class K among the N overall regions and L the number of bins of\n"
                                   "the each region.\n"
                                   "\n"
                                   "In the case of a soft clustering methods, this procedure creates a matrix of\n"
                                   "dimensions NxL' where L'=L-S+1 and S is the shifting freedom allowed during\n"
                                   "the classification. The resulting matrix contains as many rows as the\n"
                                   "starting matrix because in soft clustering, all regions (rows) are\n"
                                   "assigned to all classes\n"
                                   "\n"
                                   "To construct a final matrix M3 of dimensions NxL3 where L3 covers a given\n"
                                   "range <from>/<to>, the original matrix M1 of dimensions NxL is computed and\n"
                                   "extended into a matrix M2 NxL2 with L2>=L1. The final M3 of dimensions NxL\n"
                                   "is eventually computed, for class K, using the given posterior probabilities.\n"
                                   "A row of the final matrix M3 is the weighted average of each of the S\n"
                                   "possibles slices of the corresponding row in M2. The weights used are the\n"
                                   "probabilities with which this row was assigned to class K, for each of\n"
                                   "the S shift states, in each flip state.\n"
                                   "\n"
                                   "The original matrix M1 that was partitionned with shifting freedom S is\n"
                                   "generated using the BED, BAM and BAI files that were originally used to\n"
                                   "create it.\n"
                                   "The posterior probabilities should be a 4D matrix in binary format, with\n"
                                   "dimensions :\n"
                                   "1) number of regions\n"
                                   "2) number of classes\n"
                                   "3) number of shift states\n"
                                   "4) number of flip states\n"
                                   "The results is returned through stdout as a 2D text matrix of dimensions :\n"
                                   "1) number of regions\n"
                                   "2) number of bins in each regions, as defined by the <from>/<to> range and\n"
                                   "the <binsiz> value.\n\n" ;

    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_bed_msg      = "The path to the BED file containing the references.";
    std::string opt_bam_msg      = "The path to the BAM file containing the targets.";
    std::string opt_bai_msg      = "The path to the BAI file containing the BAM file index.";
    std::string opt_prob_msg     = "The path to the file containing the assignment probabilities\n"
                                   "(the partition)." ;
    std::string opt_from_msg     = "The most upstream position - in relative coordinate - of the regions\n"
                                   "in the original matrix." ;
    std::string opt_to_msg       = "The most downstream position - in relative coordinate - of the regions\n"
                                   "in the original matrix." ;
    std::string opt_classk_msg   = "The index (1-based) of the class of interest." ;
    std::string opt_binsize_msg  = "The size of the bins, in the original matrix." ;
    char tmp[4096] ;
    sprintf(tmp,
                                   "How the data in the BAM file should be handled when computing\n"
                                   "the number of counts in each bin.\n"
                                   "\t\"%s\" uses each position within the reads (by default)\n"
                                   "\t\"%s\" uses only the insertion site for ATAC-seq data\n"
                                   "\t\"%s\" uses each position within the fragments\n"
                                   "\t\"%s\" uses only the fragment central positions\n",
            method_read.c_str(),
            method_read_atac.c_str(),
            method_fragment.c_str(),
            method_fragment_center.c_str()) ;

     std::string opt_method_msg = tmp ;

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    std::string method(method_read) ;

    desc.add_options()
                ("help,h",  opt_help_msg.c_str())

                ("bed",     po::value<std::string>(&(this->file_bed)),  opt_bed_msg.c_str())
                ("bam",     po::value<std::string>(&(this->file_bam)),  opt_bam_msg.c_str())
                ("bai",     po::value<std::string>(&(this->file_bai)),  opt_bai_msg.c_str())
                ("prob",    po::value<std::string>(&(this->file_prob)), opt_prob_msg.c_str())

                ("from",    po::value<int>(&(this->from)),             opt_from_msg.c_str())
                ("to",      po::value<int>(&(this->to)),               opt_to_msg.c_str())
                ("binSize", po::value<int>(&(this->bin_size)),         opt_binsize_msg.c_str())
                ("k",       po::value<size_t>(&(this->class_k)),       opt_classk_msg.c_str())
                ("method",  po::value<std::string>(&(method)),         opt_method_msg.c_str()) ;

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
    else if(this->file_bam == "" and (not help))
    {   std::string msg("Error! No BAM file was given (--bam)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->file_bai == "" and (not help))
    {   std::string msg("Error! No BAM index file was given (--bai)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->file_prob == "" and (not help))
    {   std::string msg("Error! No probability (partition) file was given (--prob)!") ;
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
    else if(this->bin_size <= 0 and (not help))
    {   std::string msg("Error! bin size should be bigger than 0 (--binSize)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(method != method_read and
            method != method_read_atac and
            method != method_fragment and
            method != method_fragment_center)
    {   char msg[4096] ;
        sprintf(msg, "Error! method should be %s, %s, %s or %s (--method)",
                method_read.c_str(),
                method_read_atac.c_str(),
                method_fragment.c_str(),
                method_fragment_center.c_str()) ;
        throw std::invalid_argument(msg) ;
    }

    // set method
    if(method == method_read)
    {   this->method = CorrelationMatrixCreator::READ ; }
    else if(method == method_read_atac)
    {   this->method = CorrelationMatrixCreator::READ_ATAC ; }
    else if(method == method_fragment)
    {   this->method = CorrelationMatrixCreator::FRAGMENT ; }
    else if(method == method_fragment_center)
    {   this->method = CorrelationMatrixCreator::FRAGMENT_CENTER ; }

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
{   ClassReadDataCreatorApplication app(argn, argv) ;
    return app.run() ;
}

