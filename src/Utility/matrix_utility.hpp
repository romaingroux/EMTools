#ifndef MATRIXUTILITY_HPP
#define MATRIXUTILITY_HPP

#include <vector>

#include <Matrix2D.hpp>
#include <Matrix3D.hpp>
#include <Matrix4D.hpp>


/*!
 * \brief Filters out the given rows of a matrix.
 * \param indices the 0-based indices of the
 * rows to remove. The indices should be
 * sorted in ascending order.
 * \param m the matrix on which this operation
 * should be performed.
 * \return the matrix without the filtered
 * rows.
 * \throw std::invalid_argument if one index is
 * invalid.
 */
Matrix2D<int> filter_rows(const std::vector<size_t>& indices,
                          const Matrix2D<int>& m) ;

/*!
 * \brief Filters out the given rows (1st dimension
 * slices) of a matrix.
 * \param indices the 0-based indices of the
 * rows (1st dimension slices) to remove. The
 * indices should be sorted in ascending order.
 * \param m the matrix on which this operation
 * should be performed.
 * \return the matrix without the filtered
 * rows.
 * \throw std::invalid_argument if one index is
 * invalid.
 */
Matrix3D<double> filter_rows(const std::vector<size_t>& indices,
                             const Matrix3D<double>& m) ;


/*!
 * \brief Filters out the given rows (1st dimension
 * slices) of a matrix.
 * \param indices the 0-based indices of the
 * rows (1st dimension slices) to remove. The
 * indices should be sorted in ascending order.
 * \param m the matrix on which this operation
 * should be performed.
 * \return the matrix without the filtered
 * rows.
 * \throw std::invalid_argument if one index is
 * invalid.
 */
Matrix4D<double> filter_rows(const std::vector<size_t>& indices,
                             const Matrix4D<double>& m) ;


#endif // MATRIXUTILITY_HPP
