#include <ReadModelExtenderApplication.hpp>

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <utility>    // std::move()
#include <stdexcept>  // std::invalid_argument, std::runtime_error

#include <CorrelationMatrixCreator.hpp>
#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <ReadLayer.hpp>
#include <ReadModelComputer.hpp>

namespace po = boost::program_options ;


// the valid values for --method option
std::string method_read            = "read" ;
std::string method_read_atac       = "read_atac" ;
std::string method_fragment        = "fragment" ;
std::string method_fragment_center = "fragment_center" ;


ReadModelExtenderApplication::ReadModelExtenderApplication(int argn, char** argv)
    : file_bed(""), file_bam(""), file_bai(""), file_prob(""),
      from(0), to(0), ext(0), bin_size(0),
      method(CorrelationMatrixCreator::FRAGMENT), bckg_class(false),
      n_threads(0), runnable(false)
{   this->parseOptions(argn, argv) ; }

ReadModelExtenderApplication::~ReadModelExtenderApplication()
{}

int ReadModelExtenderApplication::run()
{   if(this->runnable)
    {   // extend limits
        int ext_right = this->ext/2 ;
        int ext_left  = this->ext - ext_right ;
        this->from -= ext_left ;
        this->to   += ext_right ;

        // create extended matrix
        CorrelationMatrixCreator mc(this->file_bed,
                                    this->file_bam,
                                    this->file_bai,
                                    this->from,
                                    this->to,
                                    this->bin_size,
                                    this->method) ;
        Matrix2D<int> data = mc.create_matrix() ;

        // compute model
        // Matrix4D<double> prob(this->file_prob) ;
        Matrix4D<double> prob ; prob.load(this->file_prob) ;
        if(prob.get_dim()[0] != data.get_nrow())
        {   char msg[4096] ;
            sprintf(msg,
                    "Error! data matrix and probability matrix have "
                    "unequal row numbers (%zu and %zu)",
                    prob.get_dim()[0],
                    data.get_nrow()) ;
            throw std::invalid_argument(msg) ;
        }
        size_t n_row   = prob.get_dim()[0] ;
        size_t n_class = prob.get_dim()[1] ;
        size_t n_shift = prob.get_dim()[2] ;
        size_t n_flip  = prob.get_dim()[3] ;

        ReadModelComputer model_cp(std::move(data),
                                   prob,
                                   this->bckg_class,
                                   this->n_threads) ;
        Matrix2D<double> model = model_cp.get_model() ;

        // compute class prob
        vector_d class_prob(n_class, 0.) ;
        double p_tot = 0. ;
        for(size_t i=0; i<n_row; i++)
        {   for(size_t j=0; j<n_class; j++)
            {   for(size_t k=0; k<n_shift; k++)
                {   for(size_t l=0; l<n_flip; l++)
                    {   class_prob[j] += prob(i,j,k,l) ;
                        p_tot         += prob(i,j,k,l) ;
                    }
                }
            }
        }
        for(auto& prob : class_prob)
        {   prob /= p_tot ; }

        // create a matrix containing the class prob in the 1st
        // column and the model in the remaining columns
        Matrix2D<double> model_final(model.get_nrow(),
                                     model.get_ncol() + 1) ;
        // 1st column contain the class prob
        for(size_t i=0; i<model_final.get_nrow(); i++)
        {   model_final(i,0) = class_prob[i] ; }

        // fill the remaining with the model parameters
        for(size_t i=0; i<model.get_nrow(); i++)
        {   for(size_t j=0; j<model.get_ncol(); j++)
            {   model_final(i,j+1) = model(i,j) ; }
        }
        std::cout << model_final << std::endl ;
        return 0 ;
    }
    else
    {   return 1 ; }
}

void ReadModelExtenderApplication::parseOptions(int argn, char** argv)
{   // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "ReadModelExtender is an application that extends a read density model of'\n"
                                   "length L' (L' = L - S + 1 where L is the number of column of the data\n"
                                   "matrix that was classified to create the models and S the shifting freedom\n"
                                   "allowed during the classification) to a new model of length L'' = L' + E\n"
                                   "(where E is the number of columns to add to the model).\n"
                                   "The extension procedure requires to extend the orignal read density matrix\n"
                                   "that was partitioned by adding E/2 columns on each side. For this, the\n"
                                   "BED, BAM and BAI files that have been used to create the original read\n"
                                   "density matrix are needed.\n"
                                   "Eventually the extended class models are generated from the extended read\n"
                                   "density matrix, using the posterior probabilites, and returned through\n"
                                   "the stdout as a 2D text matrix.\n"
                                   "The posterior probability matrix should be stored in a 4D matrix in binary\n"
                                   "format, with dimensions :\n"
                                   "1) number of regions\n"
                                   "2) number of classes\n"
                                   "3) number of shift states\n"
                                   "4) number of flip states\n"
                                   "The models are returned as a text 2D matrix of dimensions :\n"
                                   "1) number of classes\n"
                                   "2) L''+1, the 1st column contains the class probabilities\n\n" ;


    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_thread_msg   = "The number of threads dedicated to parallelize the computations,\n "
                                   "by default 0 (no parallelization)." ;
    std::string opt_bed_msg      = "The path to the BED file containing the references of the original matrix.";
    std::string opt_bam_msg      = "The path to the BAM file containing the targets of the original matrix.";
    std::string opt_bai_msg      = "The path to the BAI file containing the BAM file index.";
    std::string opt_prob_msg     = "The path to the file containing the posterior probabilities\n"
                                   "of the original matrix classification." ;
    std::string opt_from_msg     = "The upstream limit - in relative coordinate - of the region to build "
                                   "around each reference center, as it was for the original matrix." ;
    std::string opt_to_msg       = "The downstream limit - in relative coordinate - of the region to build "
                                   "around each reference center, as it was for the original matrix." ;
    std::string opt_ext_msg      = "The number of columns (E) to add to the model. Each side will be\n"
                                   "extended by E/2." ;
    std::string opt_bckg_msg     = "Whether the last class of the classification (posterior\n"
                                   "probabilities) is a background class." ;
    std::string opt_binsize_msg  = "The size of the bins, as it was for the original matrix." ;
    char tmp[4096] ;
    sprintf(tmp,
                                   "How the data in the BAM file should be handled when computing\n"
                                   "the number of counts in each bin, as it was for the original\n"
                                   "matrix.\n"
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
                ("help,h",       opt_help_msg.c_str())

                ("bed",     po::value<std::string>(&(this->file_bed)),   opt_bed_msg.c_str())
                ("bam",     po::value<std::string>(&(this->file_bam)),   opt_bam_msg.c_str())
                ("bai",     po::value<std::string>(&(this->file_bai)),   opt_bai_msg.c_str())
                ("prob,",   po::value<std::string>(&(this->file_prob)),  opt_prob_msg.c_str())

                ("bgclass", opt_bckg_msg.c_str())

                ("from",    po::value<int>(&(this->from)),               opt_from_msg.c_str())
                ("to",      po::value<int>(&(this->to)),                 opt_to_msg.c_str())
                ("ext",     po::value<int>(&(this->ext)),                opt_ext_msg.c_str())
                ("binSize", po::value<int>(&(this->bin_size)),           opt_binsize_msg.c_str())
                ("method",  po::value<std::string>(&(method)),           opt_method_msg.c_str())

                ("thread",  po::value<std::size_t>(&(this->n_threads)),  opt_thread_msg.c_str()) ;

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
    {   std::string msg("Error! No posterior probability file was given (--prob)!") ;
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
    else if(ext <= 0 and (not help))
    {   std::string msg("Error! the number of columns to add should be > 0 (--ext)!") ;
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

    // set background class
    if(vm.count("bgclass"))
    {   this->bckg_class  = true ; }

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
{   ReadModelExtenderApplication app(argn, argv) ;
    return app.run() ;
}
