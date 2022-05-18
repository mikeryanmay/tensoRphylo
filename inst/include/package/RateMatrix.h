#ifndef RMAT_H
#define RMAT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Interface/DistributionHandler.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

namespace TensorPhyloUtils {
  bool isRateMatrix(const MatrixXd& mat);
}

class RateMatrix {

  public:

    RateMatrix(size_t dim_) : dim(dim_), data( MatrixXd::Zero(dim_, dim_) ) {
    }

    RateMatrix(size_t dim_, double rate) : dim(dim_), data( rate * MatrixXd::Ones(dim_, dim_) / double(dim - 1) ) {

      // correct the diagonal value
      data.diagonal() = -rate * VectorXd::Ones(dim);

    }

    RateMatrix(const MatrixXd& mat) {

      // do some simple checking
      bool valid = TensorPhyloUtils::isRateMatrix(mat);
      if ( valid == false ) {
        stop("Error setting rate matrix. Make sure the provided matrix is a valid rate matrix.");
      }

      // get the dimensions
      dim = mat.rows();

      // set the data
      data = mat;

    }

    RateMatrix(const RateMatrix& other) {
      dim  = other.dim;
      data = other.data;
    }

    ~RateMatrix(){};

    const MatrixXd& getMatrix() const {
      return data;
    }

    const double& getRate(size_t i, size_t j) {
      return data(i,j);
    }

    void setRate(size_t i, size_t j, double val) {

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

    const size_t& getDim() const {
      return dim;
    }

    void show() {

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

  private:

    // private data
    size_t   dim;
    MatrixXd data;

};

class RateMatrixList {

  public:

    RateMatrixList(size_t i){};
    ~RateMatrixList(){};

    void addRateMatrix(RateMatrix& events) {
      matrix_list.push_back(&events);
    }

    RateMatrix& getRateMatrixOneIndex(size_t i) {
      return *(matrix_list.at(i - 1));
    }

    size_t size() const {
      return matrix_list.size();
    }

    const RateMatrix& at(size_t i) const {
      return *(matrix_list.at(i));
    }

    const std::vector<RateMatrix*>& getMatrixList() const {
      return matrix_list;
    }

    void show() {
      Rcout << "A rate matrix list with elements: " << std::endl;
      if ( matrix_list.size() > 0 ) {
        // iterate over elements
        for(std::vector<RateMatrix*>::iterator it = matrix_list.begin(); it != matrix_list.end(); ++it ) {
          (*it)->show();
        }
      }
    }

  private:

    // the vector
    std::vector<RateMatrix*> matrix_list;

};


#endif
