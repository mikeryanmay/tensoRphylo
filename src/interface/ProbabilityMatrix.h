#ifndef PMAT_H
#define PMAT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Interface/DistributionHandler.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

RCPP_EXPOSED_CLASS(ProbabilityMatrix)
class ProbabilityMatrix {

  public:

    ProbabilityMatrix(size_t dim_);
    ProbabilityMatrix(const MatrixXd& mat);
    ProbabilityMatrix(const ProbabilityMatrix& other);
    ~ProbabilityMatrix(){};

    const MatrixXd& getMatrix() const;
    const double&   getProbability(size_t i, size_t j);
    void            setProbability(size_t i, size_t j, double val);
    const size_t&   getDim() const;
    void            show();

  private:

    // private data
    size_t   dim;
    MatrixXd data;

};

RCPP_EXPOSED_CLASS(ProbabilityMatrixList)
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
