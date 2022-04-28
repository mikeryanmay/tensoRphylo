#ifndef RMAT_H
#define RMAT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include "Interface/DistributionHandler.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

// typedef std::map< std::vector<unsigned>, double > eventMap_t;

RCPP_EXPOSED_CLASS(RateMatrix)
class RateMatrix {

  public:

    RateMatrix(size_t dim_);
    RateMatrix(size_t dim_, double rate);
    RateMatrix(const RateMatrix &other);
    ~RateMatrix(){};

    const MatrixXd& getMatrix() const;
    const double&   getRate(size_t i, size_t j);
    void            setRate(size_t i, size_t j, double val);
    const size_t&   getDim() const;
    void            show();

  private:

    // private data
    size_t   dim;
    MatrixXd data;

};

RCPP_EXPOSED_CLASS(RateMatrixList)
class RateMatrixList {

  public:

    RateMatrixList(size_t i){};
    ~RateMatrixList(){};

    void addRateMatrix(RateMatrix& events) {
      matrix_list.push_back(&events);
    }

    RateMatrix getRateMatrix(size_t i) {
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
      Rcout << "A rate matrix array with elements: " << std::endl;
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
