#include <Rcpp.h>
#include <limits>
#include "RateMatrix.h"

using namespace Rcpp;
using namespace Eigen;

RateMatrix::RateMatrix(size_t dim_) : dim(dim_), data( MatrixXd::Zero(dim_, dim_) ) {

}

RateMatrix::RateMatrix(size_t dim_, double rate) : dim(dim_), data( rate * MatrixXd::Ones(dim_, dim_) / double(dim - 1) ) {

  // correct the diagonal value
  data.diagonal() = -rate * VectorXd::Ones(dim);

}

RateMatrix::RateMatrix(const MatrixXd& mat) {

  // do some simple checking
  if ( mat.rows() == mat.cols() ) {
    stop("Error copying rate matrix. Rate matrix must be square.");
  }

  // get the dimensions
  dim = mat.rows();

  // set the data
  data = mat;

}

RateMatrix::RateMatrix(const RateMatrix& other) {
  dim  = other.dim;
  data = other.data;
}

const MatrixXd& RateMatrix::getMatrix() const {
  return data;
}

const double& RateMatrix::getRate(size_t i, size_t j) {
  return data(i,j);
}

void RateMatrix::setRate(size_t i, size_t j, double val) {

  // do not allow access to the away rates
  if ( i == j ) {
    stop("Please do not attempt to access the diagonal elements.");
  }

  // do not allow index greater than dimensions
  if ( i >= dim || j >= dim ) {
    stop("Invalid index (too large).");
  }

  // make sure rates are positive
  if ( val < 0 ) {
    stop("Non-diagonal rates must be positive.");
  }

  // compute the delta
  double delta = val - data(i,j);

  // set the new value
  data(i,j) = val;

  // update the diagonal
  data(i,i) -= delta;

}


const size_t& RateMatrix::getDim() const {
  return dim;
}

void RateMatrix::show() {

  // the header
  Rcout.precision(3);
  Rcout << "A rate matrix with " << dim << " states. <" << this << ">" << std::endl;

  // loop over rows
  for(size_t i = 0; i < dim; ++i) {

    // loop over columns
    Rcout << "    [" << std::fixed << data(i,0);
    for(size_t j = 1; j < dim; ++j) {
      Rcout << ", " << data(i, j);
    } // end loop over columns
    Rcout << "]" << std::endl;

  } // end loop over rows

}











// end RateMatrix class
