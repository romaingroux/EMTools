#include <ProbToModelApplication.hpp>

#include <boost/program_options.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <utility>    // std::move()
#include <stdexcept>  // std::invalid_argument, std::runtime_error

#include <ModelComputer.hpp>
#include <ReadModelComputer.hpp>
#include <SequenceModelComputer.hpp>
#include <ConsensusSequenceModelComputer.hpp>
#include <Matrix2D.hpp>
#include <Matrix4D.hpp>
#include <matrix_utility.hpp>  // filter_rows()

namespace po = boost::program_options ;

typedef std::vector<double> vector_d ;


ProbToModelApplication::ProbToModelApplication(int argn, char** argv)
    : file_read(""), file_seq(""), file_prob(""), file_filter(""),
      bckg_class(false),
      n_threads(0), runnable(false)
{   this->parseOptions(argn, argv) ; }

ProbToModelApplication::~ProbToModelApplication()
{}

int ProbToModelApplication::run()
{   if(this->runnable)
    {
        // load data
        std::string file_data ;
        bool read_data    = false ;
        bool seq_data     = false ;
        bool consseq_data = false ;
        if(this->file_read != "")
        {   file_data = this->file_read ;
            seq_data = consseq_data = false ; read_data = true ;
        }
        else if(this->file_seq != "")
        {   file_data = this->file_seq ;
            read_data = consseq_data = false ; seq_data = true ;
        }
        else if(this->file_consseq != "")
        {   file_data = this->file_consseq ;
            read_data = seq_data = false ; consseq_data = true ;
        }
        else
        {   std::string msg("Error! Could not determine the type of the data!") ;
            throw std::runtime_error(msg) ;
        }


        // row filter
        std::vector<size_t> filter ;
        if(this->file_filter != "")
        {   // it is a column vector, easier to use the Matrix2D interface
            // to read it rather than coding a function for :)
            filter = Matrix2D<size_t>(this->file_filter).get_data() ;
            std::sort(filter.begin(), filter.end()) ;
        }

        // prob
        Matrix4D<double> prob ; prob.load(this->file_prob) ;

        // get the data model
        ModelComputer* ptr = nullptr ;
        if(read_data)
        {   Matrix2D<int> data(file_data) ;
            // filter out some rows if needed
            if(filter.size())
            {   data = filter_rows(filter, data) ; }
            this->check_dimensions(data, prob) ;
            ptr = new ReadModelComputer(std::move(data),
                                        prob,
                                        this->bckg_class,
                                        this->n_threads) ;
        }
        else if(seq_data)
        {
            Matrix2D<int> data(file_data) ;
            // filter out some rows if needed
            if(filter.size())
            {   data = filter_rows(filter, data) ; }
            this->check_dimensions(data, prob) ;
            ptr = new SequenceModelComputer(std::move(data),
                                            prob,
                                            this->bckg_class,
                                            this->n_threads) ;
        }
        else if(consseq_data)
        {   Matrix3D<double> data ;
            data.load(file_data) ;
            // filter out some rows if needed
            if(filter.size())
            {   data = filter_rows(filter, data) ; }
            this->check_dimensions(data, prob) ;
            ptr = new ConsensusSequenceModelComputer(std::move(data),
                                                     prob,
                                                     this->bckg_class,
                                                     this->n_threads) ;
        }

        Matrix2D<double> model = ptr->get_model() ;
        delete ptr ;
        ptr = nullptr ;

        // compute the class prob
        size_t n_row   = prob.get_dim()[0] ;
        size_t n_class = prob.get_dim()[1] ;
        size_t n_shift = prob.get_dim()[2] ;
        size_t n_flip  = prob.get_dim()[3] ;

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
        if(read_data)
        {   for(size_t i=0; i<model_final.get_nrow(); i++)
            {   model_final(i,0) = class_prob[i] ; }
        }
        else if(seq_data or consseq_data)
        {   size_t i_class = 0 ;
            for(size_t i=0; i<model_final.get_nrow(); i++)
            {   if((i != 0) and (i % 4 == 0))
                {   i_class++ ; }
                model_final(i,0) = class_prob[i_class] ;
            }
        }
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

void ProbToModelApplication::parseOptions(int argn, char** argv)
{   // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "ProbToModel computes the class models from\n"
                                   "the classification results (the posterior \n"
                                   "class assignment probabilities computed by\n"
                                   "EMRead, EMSequence or EMJoint and the data).\n\n"
                                   "The posterior probability should be contained in\n"
                                   "in binary format in a 4D matrix of dimensions\n"
                                   "1) number of regions\n"
                                   "2) number of classes\n"
                                   "3) number of shift states\n"
                                   "4) number offlip states)\n"
                                   "The class models is returned through stdout as a\n"
                                   "2D text matrix of dimensions :\n"
                                   "1) number of classes\n"
                                   "2) length of the models + 1, the 1st column contains\n"
                                   "the class probabilities\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_thread_msg   = "The number of threads dedicated to parallelize the computations,\n "
                                   "by default 0 (no parallelization)." ;
    std::string opt_read_msg     = "The path to the file containing the data, if the data are\n"
                                   "read counts.";
    std::string opt_seq_msg      = "The path to the file containing the data, if the data are\n"
                                   "sequences.";
    std::string opt_consseq_msg  = "The path to the file containing the data, if the data are\n"
                                   "consensus sequences.";
    std::string opt_prob_msg     = "The path to the file containing the 4D matrix with the \n"
                                   "posterior probabilities\n." ;
    std::string opt_filter_msg   = "Optional. The path to a single column text file containing the 0-based\n"
                                   "indices of rows to filter out in the data. Once the rows have been\n"
                                   "filtered out from the data, its rows order should correspond to the\n"
                                   "probabilities." ;
    std::string opt_bckg_msg     = "Whether the last class of the classification (posterior\n"
                                   "probabilities) is a background class." ;
    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    std::string seeding_tmp ;

    desc.add_options()
                ("help,h",       opt_help_msg.c_str())

                ("read,",       po::value<std::string>(&(this->file_read)),    opt_read_msg.c_str())
                ("seq,",        po::value<std::string>(&(this->file_seq)),     opt_seq_msg.c_str())
                ("consseq,",    po::value<std::string>(&(this->file_consseq)), opt_consseq_msg.c_str())
                ("prob,",       po::value<std::string>(&(this->file_prob)),    opt_prob_msg.c_str())
                ("filter",      po::value<std::string>(&(this->file_filter)),  opt_filter_msg.c_str())

                ("bgclass",     opt_bckg_msg.c_str())

                ("thread",      po::value<std::size_t>(&(this->n_threads)),    opt_thread_msg.c_str()) ;

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
    if((this->file_read    == "") and
       (this->file_seq     == "") and
       (this->file_consseq == "") and
       (not help))
    {   std::string msg("Error! No data file was given (none of --read --seq --consseq)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if((this->file_read != "") and
            (this->file_seq  != "") and
            (not help))
    {   std::string msg("Error! --read and --seq are mutually exclusive!") ;
        throw std::invalid_argument(msg) ;
    }
    else if((this->file_read    != "") and
            (this->file_consseq != "") and
            (not help))
    {   std::string msg("Error! --read and --consseq are mutually exclusive!") ;
        throw std::invalid_argument(msg) ;
    }
    else if((this->file_seq     != "") and
            (this->file_consseq != "") and
            (not help))
    {   std::string msg("Error! --seq and --consseq are mutually exclusive!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->file_prob == "" and (not help))
    {   std::string msg("Error! No posterior probabily file was given (--prob)!") ;
        throw std::invalid_argument(msg) ;
    }

    // set background class
    if(vm.count("bgclass"))
    {   this->bckg_class  = true ; }

    // help invoked, run() cannot be invoked
    if(help)
    {	std::cout << desc << std::endl ;
        this->runnable = false ;
        return ;
    }
    // everything fine, run() can be called
    else
    {   this->runnable = true ;
        return ;
    }
}

void ProbToModelApplication::check_dimensions(const Matrix2D<int>& data,
                                              const Matrix4D<double>& prob) const
{   if(data.get_nrow() != prob.get_dim()[0])
    {   char msg[4096] ;
        sprintf(msg, "Error! data and prob matrices have unequal "
                      "row numbers (%zu / %zu)!",
                data.get_nrow(), prob.get_dim()[0]) ;
        throw std::runtime_error(msg) ;
    }
    else if(data.get_ncol() < prob.get_dim()[2])
    {   char msg[4096] ;
        sprintf(msg, "Error! too many shift states for the data width!"
                     "(%zu shift states and %zu columns in data)!",
                prob.get_dim()[2], data.get_ncol()) ;
        throw std::runtime_error(msg) ;
    }
}

void ProbToModelApplication::check_dimensions(const Matrix3D<double>& data,
                                              const Matrix4D<double>& prob) const
{   if(data.get_dim()[0] != prob.get_dim()[0])
    {   char msg[4096] ;
        sprintf(msg, "Error! data and prob matrices have unequal "
                      "row numbers (%zu / %zu)!",
                data.get_dim()[0], prob.get_dim()[0]) ;
        throw std::runtime_error(msg) ;
    }
    else if(data.get_dim()[1] < prob.get_dim()[2])
    {   char msg[4096] ;
        sprintf(msg, "Error! too many shift states for the data width!"
                     "(%zu shift states and %zu data width)!",
                prob.get_dim()[2], data.get_dim()[1]) ;
        throw std::runtime_error(msg) ;
    }
    else if(data.get_dim()[2]!= 4)
    {   char msg[4096] ;
        sprintf(msg, "Error! data 3rd dimension is not 4!"
                     "(%zu)!",
                data.get_dim()[2]) ;
        throw std::runtime_error(msg) ;
    }
}


int main(int argn, char** argv)
{   ProbToModelApplication app(argn, argv) ;
    return app.run() ;
}
