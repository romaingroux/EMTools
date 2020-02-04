#ifndef MODELCOMPUTER_HPP
#define MODELCOMPUTER_HPP

#include <Matrix2D.hpp>
#include <DataLayer.hpp>

/*!
 * \brief The ModelComputer class provides the
 * bases to compute the class models, as computed
 * by one of the expectation maximization algorithm,
 * given a data matrix and a partitition.
 * The data partition must be a probability
 * matrix, as returned by one of the classes
 * deriving EMBase.
 */
class ModelComputer
{
    public:
        /*!
         * \brief Constructs an empty object.
         */
        ModelComputer() ;

        ModelComputer(const ModelComputer& other) = delete ;

        /*!
         * \brief Destructor.
         */
        virtual ~ModelComputer() ;

        /*!
         * \brief Returns the data model in a nice
         * format.
         * 1st dim: the different classes and
         * the model categories. For instance,
         * a read model with 2 classes will have
         * class 1 and class 2 over the rows.
         * A sequence model with 2 classes will
         * have class 1 A, class 1 C, class 1 G,
         * class 1 T, class 2 A, class 2 C,
         * class 2 G and class 2 T.
         * 2nd dim: the model length
         *    ___________
         *    |  class1  | /|\
         * ___|__________|_\|/ 1 (reads) or 4 (sequences)
         *    |  class2  | /|\
         *    |__________| \|/ 1 (reads) or 4 (sequences)
         *
         *    <---------->
         *     model length
         * \return the data model.
         */
        virtual Matrix2D<double> get_model() const ;

    protected:
        /*!
         * \brief The data layer containing the
         * data and their models.
         */
        DataLayer* data_layer ;
} ;

#endif // MODELCOMPUTER_HPP
