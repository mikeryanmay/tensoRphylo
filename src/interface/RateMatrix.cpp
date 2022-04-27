#include <Rcpp.h>
#include <limits>
#include "RateMatrix.h"

using namespace Rcpp;
using namespace Eigen;

RateMatrix::RateMatrix(size_t dim_) : dim(dim_), data( MatrixXd::Zero(dim_, dim_) ) {

}

RateMatrix::RateMatrix(size_t dim_, double rate) : dim(dim_), data( rate * MatrixXd::Ones(dim_, dim_) ) {

  // correct the diagonal value
  data.diagonal() = -rate * double(dim - 1) * VectorXd::Ones(dim);

}

RateMatrix::RateMatrix(const RateMatrix &other) {
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
  Rcout << "A rate matrix array with " << dim << " states. <" << this << ">" << std::endl;

  // if ( data.size() > 0 ) {
  //   // iterate over elements
  //   Rcout << "    (ancestral state -> left daughter state, right daughter state) = value" << std::endl;
  //   for(eventMap_t::iterator it = data.begin(); it != data.end(); ++it) {
  //     Rcout << "    (" << it->first[0] + 1 << " -> " << it->first[1] + 1 << ", " << it->first[2] + 1 << ") = " << it->second << std::endl;
  //   }
  // }

}











// end RateMatrix class
