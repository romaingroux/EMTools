#include <MatrixBinToTxtApplication.hpp>

#include <boost/program_options.hpp>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>


namespace po = boost::program_options ;


// valid types in the matrix handled
std::string type_int    = "int" ;
std::string type_char   = "char" ;
std::string type_double = "double" ;


MatrixBinToTxtApplication::MatrixBinToTxtApplication(int argn, char** argv)
    : file_matrix(""), type(""), ndim(0), runnable(false)
{
    // parse command line options and set the fields
    this->parseOptions(argn, argv) ;
}

int MatrixBinToTxtApplication::run()
{   if(this->runnable)
    {
        if(this->type == type_int)
        {   if(this->ndim == 2)
            {   Matrix2D<int> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
            else if(this->ndim == 3)
            {   Matrix3D<int> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
            else if(this->ndim == 4)
            {   Matrix4D<int> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
        }
        else if(this->type == type_char)
        {   if(this->ndim == 2)
            {   Matrix2D<char> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
            else if(this->ndim == 3)
            {   Matrix3D<char> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
            else if(this->ndim == 4)
            {   Matrix4D<char> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
        }
        else if(this->type == type_double)
        {   if(this->ndim == 2)
            {   Matrix2D<double> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
            else if(this->ndim == 3)
            {   Matrix3D<double> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
            else if(this->ndim == 4)
            {   Matrix4D<double> m ; m.load(this->file_matrix) ; std::cout << m << std::endl ; }
        }
        return EXIT_SUCCESS ;
    }
    else
    {   return EXIT_FAILURE; }
}

void MatrixBinToTxtApplication::parseOptions(int argn, char** argv)
{
    // no option to parse
    if(argv == nullptr)
    {   std::string message = "no options to parse!" ;
        throw std::invalid_argument(message) ;
    }

    // help messages
    std::string desc_msg =         "\n"
                                   "MatrixBinToTxt reads a Matrix that has been dumped in file in\n"
                                   "binary format (using Matrix.save(const std::string&) method)and "
                                   "displays it in text format through stdout.\n\n" ;
    std::string opt_help_msg     = "Produces this help message." ;
    std::string opt_file_mat_msg = "The path to the matrix file." ;
    std::string opt_ndim_msg     = "The number of dimensions that the matrix has." ;
    char opt_type_msg[4096] ;
    sprintf(opt_type_msg,          "The type of the data stored in the matrix. It should be %s, %s or\n"
                                   "%s.", type_int.c_str(), type_char.c_str(), type_double.c_str());

    // option parser
    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    desc.add_options()
                ("help,h",   opt_help_msg.c_str())

                ("file",     po::value<std::string>(&(this->file_matrix)), opt_file_mat_msg.c_str())
                ("type",     po::value<std::string>(&(this->type)),        opt_type_msg)
                ("ndim",     po::value<size_t>(&(this->ndim)),             opt_ndim_msg.c_str()) ;

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
    {   std::string msg("Error! No input file was given (--file)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->type == "" and
            (not help))
    {   std::string msg("Error! No data type was not given (--type)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->ndim == 0 and
            (not help))
    {   std::string msg("Error! The matrix dimensionality was not given (--ndim)!") ;
        throw std::invalid_argument(msg) ;
    }
    else if((this->type != type_int)    and
            (this->type != type_char)   and
            (this->type != type_double) and
            (not help))
    {   char msg[4096] ;
        sprintf(msg, "Error! Unsupported data type %s! It should be %s, %s or %s (--type)",
                this->type.c_str(),
                type_int.c_str(), type_char.c_str(), type_double.c_str()) ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->ndim == 1)
    {   std::string msg("Error! The matrix dimensionality cannot be 1!") ;
        throw std::invalid_argument(msg) ;
    }
    else if(this->ndim > 4)
    {   char msg[4096] ;
        sprintf(msg, "Error! The matrix dimensionality cannot be bigger than 4 (%zu)!",
                this->ndim) ;
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
{   MatrixBinToTxtApplication app(argn, argv) ;
    return app.run() ;
}
