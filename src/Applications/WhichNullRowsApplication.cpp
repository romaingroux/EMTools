#include <WhichNullRowsApplication.hpp>

#include <vector>

#include <boost/program_options.hpp>

#include <Matrix2D.hpp>


namespace po = boost::program_options ;



WhichNullRowsApplication::WhichNullRowsApplication(int argn, char** argv)
    :file_matrix(""), runnable("")
{   this->parseOptions(argn, argv) ; }

int WhichNullRowsApplication::run()
{
    if(this->runnable)
    {   Matrix2D<int> m(this->file_matrix) ;

        std::vector<size_t> null_rows ;

        for(size_t i=0; i<m.get_nrow(); i++)
        {   bool null_row = true ;
            for(size_t j=0; j<m.get_ncol(); j++)
            {   // this row has signal
                if(m(i,j) != 0)
                {   null_row = false ;
                    break ;
                }
            }
            if(null_row)
            {   null_rows.push_back(i) ; }
        }

        // clean memory
        m = Matrix2D<int>() ;

        // put the vector into a single column matrix
        // to use the << operator
        Matrix2D<size_t> rows(null_rows.size(), 1) ;
        for(size_t i=0; i< null_rows.size(); i++)
        {   rows(i,0) = null_rows[i] ; }

        std::cout << rows << std::endl ;

        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE; }
}

void WhichNullRowsApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "WhichNullRow reads a 2D matrix file in text format "
                                   "containing read densities and returns the index (0-based)\n"
                                   "of the rows which contain no read. The dimensions\n"
                                   "of the matrix should be :\n"
                                   "1) number of regions\n"
                                   "2) number of bins per region\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_file_mat_msg = "The path to the 2D matrix file." ;

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    desc.add_options()
                ("help,h",   opt_help_msg.c_str())

                ("mat",     po::value<std::string>(&(this->file_matrix)), opt_file_mat_msg.c_str()) ;

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
    if(this->file_matrix == "" and
       (not help))
    {   std::string msg("Error! No input file was given (--mat)!") ;
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
{   WhichNullRowsApplication app(argn, argv) ;
    return app.run() ;
}
