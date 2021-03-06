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

inline stdVectorXd ScalarToStdVec(const double& scalar, size_t dim) {
  return stdVectorXd(dim, scalar);
}

inline stdMatrixXd ScalarToStdMat(const double& scalar, size_t dim) {
  return stdMatrixXd(1, stdVectorXd(dim, scalar));
}

inline stdMatrixXd VectorToStdMat(const VectorXd& vec, size_t dim) {

  // copy the data
  stdVectorXd row(vec.data(), vec.data() + vec.size());
  stdMatrixXd res(dim, row);

  return res;

}

inline stdVectorXd EigenToStd(const VectorXd& eig_vec) {

  // create the vector
  stdVectorXd res(eig_vec.data(), eig_vec.data() + eig_vec.size());

  return res;

}

inline stdMatrixXd EigenToStd(const MatrixXd& eig_mat) {

  // create the matrix
  stdMatrixXd mat;
  for(size_t r = 0; r < eig_mat.rows(); ++r) {
    stdVectorXd row;
    for(size_t c = 0; c < eig_mat.cols(); ++c) {
      row.push_back( eig_mat(r, c) );
    }
    mat.push_back(row);
  }

  return mat;

}

inline stdMatrixXd TransProbToStd(const ProbabilityMatrix& prob_mat) {
  // get the internal matrix and convert it
  return EigenToStd(prob_mat.getMatrix());
}

inline stdMatrixXd RateMatToStd(const RateMatrix& rate_mat) {
  // get the internal matrix and convert it
  return EigenToStd(rate_mat.getMatrix());
}


inline std::vector<stdMatrixXd> EigenToStd(const std::vector<MatrixXd>& eig_mats) {

  // make the empty vector
  std::vector<stdMatrixXd> res(eig_mats.size());

  // emplace each value
  for(size_t i = 0; i < res.size(); ++i) {
    res.at(i) = EigenToStd(eig_mats.at(i));
  }

  return res;

}

inline std::vector<stdMatrixXd> TransProbMatrixToStd(const ProbabilityMatrixList& prob_mats) {

  // make the empty vector
  std::vector<stdMatrixXd> res(prob_mats.size());

  // emplace each value
  for(size_t i = 0; i < res.size(); ++i) {
    res.at(i) = TransProbToStd(prob_mats.at(i));
  }

  return res;

}

inline std::vector<stdMatrixXd> RateMatrixToStd(const RateMatrixList& rate_mats) {

  // make the empty vector
  std::vector<stdMatrixXd> res(rate_mats.size());

  // emplace each value
  for(size_t i = 0; i < res.size(); ++i) {
    res.at(i) = RateMatToStd(rate_mats.at(i));
  }

  return res;

}

inline VectorXd StdToEigen(const stdVectorXd& std_vec) {

  // create the vector
  VectorXd res(std_vec.size());
  for(size_t i = 0; i < res.size(); ++i) {
    res(i) = std_vec.at(i);
  }

  return res;

}

inline MatrixXd StdToEigen(const stdMatrixXd& std_mat) {

  // return an empty matrix
  if ( std_mat.size() == 0 ) {
    MatrixXd mat(0,0);
    return mat;
  }

  // create the matrix
  MatrixXd mat( std_mat.size(), std_mat[0].size() );
  for(size_t r = 0; r < mat.rows(); ++r) {
    for(size_t c = 0; c < mat.cols(); ++c) {
      mat(r,c) = std_mat.at(r).at(c);
    }
  }

  return mat;

}

////////////////
// validators //
////////////////

/////////////////////////
// check close to zero //
/////////////////////////

inline bool isCloseToZero(const double& x, double tol = 1e-10) {

  // check close to zero
  return std::abs(x) < std::numeric_limits<double>::epsilon();

}

/////////////////////////
// validate dimensions //
/////////////////////////

// make sure length is correct
inline bool hasDimensions(const stdVectorXd& vec, size_t size) {
  return vec.size() == size;
}

// make sure length is correct
inline bool hasDimensions(const VectorXd& vec, size_t size) {
  return vec.size() == size;
}

// this is a little complicated because the rows of this matrix
// need not be the same length!
inline bool hasDimensions(const stdMatrixXd& mat, size_t nrow, size_t ncol) {

  // check number of rows
  if ( mat.size() != nrow ) {
    return false;
  }

  // otherwise, we have to check each row to make
  // sure it has the right number of columns
  for(size_t r = 0; r < mat.size(); ++r) {

    // check the row for the right number of columns
    if ( mat.at(r).size() != ncol ) {
      return false;
    }

  } // end loop over rows

  // if we made it this far, we passed
  return true;

}

// make sure number of rows and columns is correct
inline bool hasDimensions(const MatrixXd& mat, size_t nrow, size_t ncol) {

  // this is easy
  return mat.rows() == nrow && mat.cols() == ncol;

}

// make sure the cube has the right number of rows, columns, and slices
inline bool hasDimensions(const std::vector<stdMatrixXd>& x, size_t nrow, size_t ncol, size_t nslice) {

  // check number of slices
  if ( x.size() != nslice ) {
    return false;
  }

  // now loop over each slice...
  for(size_t s = 0; s < x.size(); ++s) {

    // get the slice
    const stdMatrixXd& this_mat = x.at(s);

    // check the number of rows
    if ( this_mat.size() != nrow ) {
      return false;
    }

    // now loop over rows
    for(size_t r = 0; r < this_mat.size(); ++r) {
      if ( this_mat.at(r).size() != ncol ) {
        return false;
      }
    } // end loop over rows

  } // end loop over slices

  // if we made it this far, we passed
  return true;

}

// make sure number of rows and columns is correct
inline bool hasDimensions(const RateMatrix& mat, size_t nrow, size_t ncol) {
  return hasDimensions(mat.getMatrix(), nrow, ncol);
}

// make sure number of rows and columns is correct
inline bool hasDimensions(const ProbabilityMatrix& mat, size_t nrow, size_t ncol) {
  return hasDimensions(mat.getMatrix(), nrow, ncol);
}

// make sure the events have the right number of dimensions
inline bool hasDimensions(const CladoEvents& map, size_t dim) {

  // check that the dimensions are right
  return map.getDim() == dim;

}

//////////////////////////////
// validate positive values //
//////////////////////////////

// ensure that the value is greater than 0
inline bool isStrictlyNonNegative(const double& x) {
  return !(x < 0.0);
}

// ensure that all values are greater than 0
inline bool isStrictlyNonNegative(const stdVectorXd& x) {

  for(size_t i = 0; i < x.size(); ++i) {
    if ( x.at(i) < 0.0 ) {
      return false;
    }
  }

  return true;

}

// ensure that all values are greater than 0
inline bool isStrictlyNonNegative(const stdMatrixXd& x) {

  // loop over each row
  for(size_t r = 0; r < x.size(); ++r) {

    // check if the row passes
    bool row_pass = isStrictlyNonNegative(x.at(r));
    if ( row_pass == false ) {
      return false;
    }

  }

  // otherwise pass
  return true;

}

// ensure that all values are greater than 0
inline bool isStrictlyNonNegative(const VectorXd& x) {

  // loop over values
  for(size_t i = 0; i < x.size(); ++i) {
    if ( x(i) < 0.0 ) {
      return false;
    }
  }

  return true;

}

// ensure that all values are greater than 0
inline bool isStrictlyNonNegative(const MatrixXd& x) {

  // loop over rows
  for(size_t r = 0; r < x.rows(); ++r) {

    // loop over columns
    for(size_t c = 0; c < x.cols(); ++c) {
      if ( x(r,c) < 0.0 ) {
        return false;
      }
    }

  }

  return true;

}

////////////////////////
// validate simplexes //
////////////////////////

// ensure that all values are probabilities, and
// that the vector sums to one
inline bool isSimplex(const stdVectorXd& x) {

  // initialize accumulator
  double sum = 0.0;

  // loop over elements
  for(size_t i = 0; i < x.size(); ++i) {

    // if the element is not a probability, return false
    if ( x.at(i) < 0.0 || x.at(i) > 1.0 ) {
      return false;
    }

    // otherwise, increment the sum and continue
    sum += x.at(i);

  }

  // make sure sum is very close to 1
  if ( isCloseToZero(sum - 1.0) == false ) {
    return false;
  }

  // otherwise, it worked
  return true;

}


// ensure that all values are probabilities, and
// that the vector sums to one
inline bool isSimplex(const VectorXd& x) {

  // initialize accumulator
  double sum = 0.0;

  // loop over elements
  for(size_t i = 0; i < x.size(); ++i) {

    // if the element is not a probability, return false
    if ( x(i) < 0.0 || x(i) > 1.0 ) {
      return false;
    }

    // otherwise, increment the sum and continue
    sum += x(i);

  }

  // make sure sum is very close to 1
  if ( isCloseToZero(sum - 1.0) == false ) {
    return false;
  }

  // otherwise, it worked
  return true;

}

////////////////////////////
// validate probabilities //
////////////////////////////

inline bool isProbability(const double& x) {
  return !(x < 0.0 || x > 1.0);
}

inline bool isProbability(const stdVectorXd& x) {

  // loop over each element
  for(size_t i = 0; i < x.size(); ++i) {
    if ( x.at(i) < 0.0 || x.at(i) > 1.0 ) {
      return false;
    }
  }

  // otherwise, we passed
  return true;

}

inline bool isProbability(const VectorXd& x) {

  // loop over each element
  for(size_t i = 0; i < x.size(); ++i) {
    if ( x(i) < 0.0 || x(i) > 1.0 ) {
      return false;
    }
  }

  // otherwise, we passed
  return true;

}

inline bool isProbability(const stdMatrixXd& x) {

  // convert to an eigen matrix then use appropriate
  // function
  return isProbability(StdToEigen(x));

}

inline bool isProbability(const MatrixXd& x) {

  // loop over rows
  for(size_t r = 0; r < x.rows(); ++r) {

    // loop over columns
    for(size_t c = 0; c < x.cols(); ++c) {

      if ( x(r,c) < 0.0 || x(r,c) > 1.0 ) {
        return false;
      }

    } // end loop over columns

  } // end loop over rows

  // otherwise, we passed
  return true;

}

inline bool isProbability(const ProbabilityMatrix& x) {
  // just get the internal matrix and check that
  return isProbability(x.getMatrix());
}

inline bool isProbability(const CladoEvents& x) {

  // make sure every element of the event map is a probability
  const eventMap_t map = x.getEvents();
  for( eventMap_t::const_iterator it = map.begin(); it != map.end(); ++it) {

    // get check the value is a probability
    if ( isProbability(it->second) == false ) {
      return false;
    }

  }

  // otherwise, we passed
  return true;

}

////////////////////////////
// validate rate matrices //
////////////////////////////

// ensure the diagonal elements are equal to the sum of the other
// values in the corresponding row
inline bool isRateMatrix(const MatrixXd& x) {

  // make sure the dimensions are correct
  if ( x.rows() != x.cols() ) {
    return false;
  }

  // get the diagonal elements
  const VectorXd& diag = x.diagonal();

  // make sure diagonal is strictly negative
  if (isStrictlyNonNegative(diag) == true) {
    return false;
  }

  // loop over rows
  for(size_t r = 0; r < x.rows(); ++r) {

    // loop over columns
    double row_sum = 0.0;
    for(size_t c = 0; c < x.cols(); ++c) {
      if ( r != c ) {
        row_sum += x(r,c);
      }
    } // end loop over columns

    // check if this row is good
    if ( isCloseToZero(diag(r) + row_sum) == false ) {
      return false;
    }

  } // end loop over rows

  // otherwise, we passed
  return true;

}

// ensure the diagonal elements are equal to the sum of the other
// values in the corresponding row
inline bool isRateMatrix(const stdMatrixXd& x) {

  // just convert to an eigen object and call
  // appropriate function
  return isRateMatrix(StdToEigen(x));

}

// ensure the diagonal elements are equal to the sum of the other
// values in the corresponding row
inline bool isRateMatrix(const RateMatrix& x) {

  return isRateMatrix( x.getMatrix() );

}

//////////////////////////////////////////////
// validate transition probability matrices //
//////////////////////////////////////////////

inline bool isTransitionProbabilityMatrix(const MatrixXd& x) {

  // make sure the dimensions are correct
  if ( x.rows() != x.cols() ) {
    return false;
  }

  // check that the rows sum to 1
  VectorXd s = x.rowwise().sum();
  for(size_t i = 0; i < x.rows(); ++i) {
    if ( isCloseToZero(s(i) - 1.0) == false ) {
      return false;
    }
  }

  // otherwise, we passed
  return true;

}

inline bool isTransitionProbabilityMatrix(const stdMatrixXd& x) {

  // just convert to an eigen object and call
  // appropriate function
  return isTransitionProbabilityMatrix(StdToEigen(x));

}

inline bool isTransitionProbabilityMatrix(const ProbabilityMatrix& x) {

  // get the eigen matrix and dispatch appropriate method
  return isTransitionProbabilityMatrix(x.getMatrix());

}


//////////////////////////////////
// validate cladogenetic events //
//////////////////////////////////

inline bool isCladogeneticProbabilityMap(const CladoEvents& x) {

  // get the event map
  size_t dim = x.getDim();
  const eventMap_t& map = x.getEvents();

  // make sure the parental events sum to 1
  std::vector<double> sums(dim);
  for(eventMap_t::const_iterator it = map.begin(); it != map.end(); ++it) {

    // get the ancestor index
    unsigned idx = it->first[0];

    // check that the value is between 0 and 1
    if (isProbability( it->second ) == false) {
      return false;
    }

    // increment the sum
    sums.at(idx) += it->second;

  }

  // check the sums
  for(size_t i = 0; i < dim; ++i) {
    if ( isCloseToZero(sums.at(i) - 1.0) == false ) {
      return false;
    }
  }

  // if we got this far, we passed
  return true;

}


}

#endif
