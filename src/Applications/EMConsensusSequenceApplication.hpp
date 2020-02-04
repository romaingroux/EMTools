#ifndef EMCONSENSUSSEQUENCEAPPLICATION_HPP
#define EMCONSENSUSSEQUENCEAPPLICATION_HPP

#include <ApplicationInterface.hpp>

#include <string>
#include <vector>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>

/*!
 * \brief The EMConsensusSequenceApplication class is a wrapper around an
 * EMConsensusSequence instance creating an autonomous application to classify
 * DNA consensus sequences by directly passing all the options and parameters from
 * the command line.
 * A DNA consensus sequence represents a DNA sequences as probability matrix which
 * contains the probability of each base at each position. The consensus sequence
 * matrix is a 3D matrix in binary format which dimensions are :
 * 1) the number of sequences
 * 2) the length of the sequences
 * 3) 4 for A, C, G, T.
 * The probabilites for a given position (2nd dimension) in a given consensus
 * sequence (1st dimension) must sum to 1.
 * The assignment probabilities are written in binary format as a 4D matrix
 * which dimensions are :
 * 1) the number of sequences
 * 2) the number of classes
 * 3) the number of shift states
 * 4) the number of flip states
 */
class EMConsensusSequenceApplication: public ApplicationInterface
{
    public:
        EMConsensusSequenceApplication() = delete ;
        EMConsensusSequenceApplication(const EMConsensusSequenceApplication& app) = delete ;
        /*!
         * \brief Constructs an object from the command line
         * options.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         */
        EMConsensusSequenceApplication(int argn, char** argv) ;

        /*!
         * \brief Runs the application. The data are classified
         * using the given settings and the posterior probability
         * matrix is returned through the stdout.
         * The matrix is a 4D matrix with dimensions :
         * regions, class, shift flip.
         * \return an exit code  EXIT_SUCCESS or EXIT_FAILURE
         * to return to the OS.
         */
        virtual int run() override ;

    private:
        /*!
         * \brief Parses the program command line options and
         * sets the object field accordingly.
         * If the help option is detected, the "runnable"
         * field is set to false and subsequent calls to
         * run() will produce nothing.
         * \param argn the number of options passed to the
         * main() function.
         * \param argv the vector of options passed to the
         * main() function.
         * \throw std::invalid_argument if an error is found
         * in the program options.
         */
        void parseOptions(int argn, char** argv) ;

        /*!
         * \brief Initialise the class models from the
         * most significantly enriched kmers in the data.
         * The kmers are turned into probability matrices
         * using the following rules :
         * 1) the probability of a base matching a base at
         * a given position of the kmer is set to 0.7, all
         * other to 0.1
         * 2) if the model length is bigger than 4, two
         * N's are added on each side of the central kmer
         * core (prob are 0.25, 0,25, 0.25, 0.25).
         * \param l_model the length of the model.
         * \param data the sequences in integer format
         * (A:0, C:1, G:2, T:3).
         * Its dimensions are the following :
         * 1st the number of sequences
         * 2nd the length of the sequences.
         * 3rd 4 for A, C, G, T.
         * \return The models initialised.
         * Its dimensions are the following :
         * 1st the number of classes
         * 2nd the length of the model
         * 3rd 4 for A, C, G, T
         */
        Matrix3D<double> init_model_kmer(size_t l_model,
                                         const Matrix3D<double>& data) const ;

        /*!
         * \brief the paths to the file containing the
         * consensus sequence data.
         */
        std::string file_consseq ;

        /*!
         * \brief the path to the file containing the
         * (0-based) index of the rows to filter out
         * from the data.
         */
        std::string file_filter ;

        /*!
         * \brief the path to the file in which the probability
         * matrix will be saved.
         */
        std::string file_out ;

        /*!
         * \brief the number of classes to partition the data into.
         */
        size_t n_class ;
        /*!
         * \brief the number of iterations allowed.
         */
        size_t n_iter ;
        /*!
         * \brief the shifting freedom.
         */
        size_t n_shift ;
        /*!
         * \brief whether flipping freedom is allowed.
         */
        bool flip ;
        /*!
         * \brief whether a constant class to model the
         * sequence background should be added. This
         * class has the sequence background probabilities
         * at each position.
         */
        bool bckg_class ;
        /*!
         * \brief the number of threads.
         */
        size_t n_threads ;
        /*!
         * \brief a seed to initialise the random number generator.
         */
        std::string seed ;
        /*!
         * \brief a flag indicating whether the core of run() can be
         * run or not.
         */
        bool runnable ;
} ;

#endif // EMCONSENSUSSEQUENCEAPPLICATION_HPP
