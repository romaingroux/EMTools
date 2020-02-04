
#include <CorrelationMatrixCreatorApplication.hpp>
#include <CorrelationMatrixCreator.hpp>

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
    : file_bed(""), file_bam(""), file_bai(""), from(0), to(0), bin_size(0),
      method(CorrelationMatrixCreator::FRAGMENT), runnable(true)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int CorrelationMatrixCreatorApplication::run()
{   if(this->runnable)
    {   CorrelationMatrixCreator mc(this->file_bed,
                                    this->file_bam,
                                    this->file_bai,
                                    this->from,
                                    this->to,
                                    this->bin_size,
                                    this->method) ;

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
                                   "CorrelationMatrixCreator is an application that creates a\n"
                                   "count matrix from a BED file and a BAM file and returnes it\n"
                                   "through stdout.\n"
                                   "The center of each region is computed as the center of the\n"
                                   "BED regions. Then a set of equally sized, non-overlapping\n"
                                   "bins, centered on the region center and covering the interval\n"
                                   "[from,to] is build for each region. Then, each bin is assigned\n"
                                   "the number of read/fragment positions (targets) present in\n"
                                   "the BAM file that are mapped at that position.\n"
                                   "The read/fragment counts are then computed, for each bin,\n"
                                   "from the BAM file.\n"
                                   "The matrix is a 2D matrix which dimensions are :\n"
                                   "1) number of regions\n"
                                   "2) length of region (to - from + 1) / bin_size\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_bed_msg      = "The path to the BED file containing the references.";
    std::string opt_bam_msg      = "The path to the BAM file containing the targets.";
    std::string opt_bai_msg      = "The path to the BAI file containing the BAM file index.";
    std::string opt_from_msg     = "The upstream limit - in relative coordinate - of the region to build\n"
                                   "around each reference center." ;
    std::string opt_to_msg       = "The downstream limit - in relative coordinate - of the region to build\n"
                                   "around each reference center." ;
    std::string opt_binsize_msg  = "The size of the bins." ;
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

                ("bed",     po::value<std::string>(&(this->file_bed)), opt_bed_msg.c_str())
                ("bam",     po::value<std::string>(&(this->file_bam)), opt_bam_msg.c_str())
                ("bai",     po::value<std::string>(&(this->file_bai)), opt_bai_msg.c_str())

                ("from",    po::value<int>(&(this->from)),             opt_from_msg.c_str())
                ("to",      po::value<int>(&(this->to)),               opt_to_msg.c_str())
                ("binSize", po::value<int>(&(this->bin_size)),         opt_binsize_msg.c_str())
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
    else if(this->from == 0 and this->to == 0 and (not help))
    {   std::string msg("Error! No range given (--from and --to)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->from >= this->to and (not help))
    {   std::string msg("Error! from shoud be smaller than to (--from and --to)!") ;
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
{   CorrelationMatrixCreatorApplication app(argn, argv) ;
    return app.run() ;
}

