#ifndef TPUTILS_H
#define TPUTILS_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Interface/DistributionHandlerImpl.h"
#include "CladoEvents.h"
#include "RateMatrix.h"
#include "ProbabilityMatrix.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

namespace TensorPhyloUtils {

  ////////////////
  // converters //
  ////////////////

  stdVectorXd ScalarToStdVec(const double& scalar, size_t dim);
  stdMatrixXd ScalarToStdMat(const double& scalar, size_t dim);
  stdMatrixXd VectorToStdMat(const VectorXd& vec, size_t dim);

  stdVectorXd EigenToStd(const VectorXd& eig_vec);
  stdMatrixXd EigenToStd(const MatrixXd& eig_mat);
  stdMatrixXd TransProbToStd(const ProbabilityMatrix& prob_mat);
  stdMatrixXd RateMatToStd(const RateMatrix& rate_mat);

  std::vector<stdMatrixXd> EigenToStd(const std::vector<MatrixXd>& eig_mats);
  std::vector<stdMatrixXd> TransProbMatrixToStd(const ProbabilityMatrixList& prob_mats);
  std::vector<stdMatrixXd> RateMatrixToStd(const RateMatrixList& rate_mats);

  VectorXd StdToEigen(const stdVectorXd& std_vec);
  MatrixXd StdToEigen(const stdMatrixXd& std_mat);

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
  bool hasDimensions(const RateMatrix& mat, size_t nrow, size_t ncol);
  bool hasDimensions(const ProbabilityMatrix& mat, size_t nrow, size_t ncol);
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
  bool isProbability(const ProbabilityMatrix& x);
  bool isProbability(const CladoEvents& x);

  // ensure the diagonals are equal to the (non-diagonal) rowsums
  bool isRateMatrix(const stdMatrixXd& x);
  bool isRateMatrix(const MatrixXd& x);
  bool isRateMatrix(const RateMatrix& x);

  // ensure the rows sum to 1
  bool isTransitionProbabilityMatrix(const stdMatrixXd& x);
  bool isTransitionProbabilityMatrix(const MatrixXd& x);
  bool isTransitionProbabilityMatrix(const ProbabilityMatrix& x);

  // ensure parent probabilities sum to 1
  bool isCladogeneticProbabilityMap(const CladoEvents& x);

}

#endif
