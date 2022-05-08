#ifndef PMAT_H
#define PMAT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Interface/DistributionHandler.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

class ProbabilityMatrix {

  public:

    ProbabilityMatrix(size_t dim_) : dim(dim_), data( MatrixXd::Zero(dim_, dim_) ) {

      // correct the diagonal value
      data.diagonal() = VectorXd::Ones(dim);

    }

    ProbabilityMatrix(const MatrixXd& mat)  {

      // do some simple checking
      if ( mat.rows() == mat.cols() ) {
        stop("Error copying probability matrix. Probability matrix must be square.");
      }

      // get the dimensions
      dim = mat.rows();

      // set the data
      data = mat;

    }

    ProbabilityMatrix(const ProbabilityMatrix& other) {
      dim  = other.dim;
      data = other.data;
    }

    ~ProbabilityMatrix(){};

    const MatrixXd& getMatrix() const {
      return data;
    }

    const double& getProbability(size_t i, size_t j)  {
      return data(i,j);
    }

    void setProbability(size_t i, size_t j, double val)  {

      // do not allow access to the away Probabilitys
      if ( i == j ) {
        stop("Please do not attempt to access the diagonal elements.");
      }

      // compute the delta
      double delta = val - data(i,j);

      // check for valid probability
      if ( val > 1.0 || val < 0.0 || (data(i,i) - delta) > 1.0 || (data(i,i) - delta) < 0.0 ) {
        stop("Elements of probability matrix must be between 0 and 1 (inclusive), and rows must sum to 1.");
      }

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
      Rcout << "A Probability matrix with " << dim << " states. <" << this << ">" << std::endl;

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

class ProbabilityMatrixList {

  public:

    ProbabilityMatrixList(size_t i){};
    ~ProbabilityMatrixList(){};

    void addProbabilityMatrix(ProbabilityMatrix& events) {
      matrix_list.push_back(&events);
    }

    ProbabilityMatrix getProbabilityMatrix(size_t i) {
      return *(matrix_list.at(i - 1));
    }

    size_t size() const {
      return matrix_list.size();
    }

    const ProbabilityMatrix& at(size_t i) const {
      return *(matrix_list.at(i));
    }

    const std::vector<ProbabilityMatrix*>& getMatrixList() const {
      return matrix_list;
    }

    void show() {
      Rcout << "A Probability matrix list with elements: " << std::endl;
      if ( matrix_list.size() > 0 ) {
        // iteProbability over elements
        for(std::vector<ProbabilityMatrix*>::iterator it = matrix_list.begin(); it != matrix_list.end(); ++it ) {
          (*it)->show();
        }
      }
    }

  private:

    // the vector
    std::vector<ProbabilityMatrix*> matrix_list;

};


#endif
