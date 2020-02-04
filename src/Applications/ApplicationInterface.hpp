
/*!
 * \brief The ApplicationInterface class provides an interface for
 * classes that allow to create autonomous application.
 * Running the application then only requires to run the run()
 * method that is declared in this class.
 */
class ApplicationInterface
{
   public:
        /*!
         * \brief Destructor.
         */
        virtual ~ApplicationInterface() ;
        /*!
         * \brief Runs the application, with all its
         * functionalities.
         * \return the exit code to return, for instance,
         * to the OS.
         */
        virtual int run() = 0 ;
} ;
