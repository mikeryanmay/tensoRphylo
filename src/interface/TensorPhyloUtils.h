#ifndef TPUTILS_H
#define TPUTILS_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "Interface/DistributionHandlerImpl.h"
#include "CladoEvents.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

namespace TensorPhyloUtils {

  ////////////////
  // converters //
  ////////////////

  stdVectorXd EigenToStd(const VectorXd& eig_vec);
  VectorXd    StdToEigen(const stdVectorXd& std_vec);

  stdMatrixXd EigenToStd(const MatrixXd& eig_mat);
  MatrixXd    StdToEigen(const stdMatrixXd& std_mat);

  std::vector<stdMatrixXd> ArmaToStd(const arma::cube& arma_cube);
  arma::cube               StdToArma(const std::vector<stdMatrixXd>& std_cube);

  ////////////////
  // validators //
  ////////////////

  // check close to zero
  bool isCloseToZero(const double& x, double tol = 1e-10);

  // ensure dimensionality
  bool hasDimensions(const stdVectorXd& vec, size_t size);
  bool hasDimensions(const VectorXd& vec, size_t size);
  bool hasDimensions(const stdMatrixXd& mat, size_t nrow, size_t ncol);
  bool hasDimensions(const MatrixXd& mat, size_t nrow, size_t ncol);
  bool hasDimensions(const std::vector<stdMatrixXd>& x, size_t nrow, size_t ncol, size_t nslice);
  bool hasDimensions(const arma::cube& x, size_t nrow, size_t ncol, size_t nslice);
  bool hasDimensions(const CladoEvents& map, size_t dim);

  // ensure that all values are >= 0
  bool isStrictlyNonNegative(const double& x);
  bool isStrictlyNonNegative(const stdVectorXd& x);
  bool isStrictlyNonNegative(const VectorXd& x);
  bool isStrictlyNonNegative(const stdMatrixXd& x);
  bool isStrictlyNonNegative(const MatrixXd& x);

  // ensure that values sum to 1
  bool isSimplex(const stdVectorXd& x);
  bool isSimplex(const VectorXd& x);

  // ensure all probabilities are at least 0
  bool isProbability(const double& x);
  bool isProbability(const stdVectorXd& x);
  bool isProbability(const VectorXd& x);
  bool isProbability(const stdMatrixXd& x);
  bool isProbability(const MatrixXd& x);
  bool isProbability(const CladoEvents& x);

  // ensure the diagonals are equal to the (non-diagonal) rowsums
  bool isRateMatrix(const stdMatrixXd& x);
  bool isRateMatrix(const MatrixXd& x);
  bool isRateMatrix(const arma::mat& x);

  // ensure parent probabilities sum to 1
  bool isCladogeneticProbabilityMap(const CladoEvents& x);

}

#endif
