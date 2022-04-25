#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "TensorPhyloExternal.h"

using namespace Rcpp;
// using namespace Eigen;
using namespace TensorPhylo::Interface;

// constructors/destructors //
TensorPhyloExternal::TensorPhyloExternal() {
  internal = DistributionHandlerImpl::create();
};

TensorPhyloExternal::TensorPhyloExternal(size_t dim_) :
  dim(dim_),
  root_frequency( stdVectorXd(dim_, 1.0 / (double)dim_) ),
  lambdas( stdMatrixXd(1, stdVectorXd(dim_, 1.0)) ),
  lambda_times( stdVectorXd() ),
  mus( stdMatrixXd(1, stdVectorXd(dim_, 0.0)) ),
  mu_times( stdVectorXd() ),
  phis( stdMatrixXd(1, stdVectorXd(dim_, 0.0)) ),
  phi_times( stdVectorXd() ),
  deltas( stdMatrixXd(1, stdVectorXd(dim_, 0.0)) ),
  delta_times( stdVectorXd() ),
  etas( std::vector<stdMatrixXd>(1, stdMatrixXd(dim_, stdVectorXd(dim_, 0.0)) ) ),
  eta_times( stdVectorXd() )
{
  // create the internal
  internal = DistributionHandlerImpl::create();

  // set the internal defaults
  internal->setRootPrior(root_frequency);
  internal->setLambda(lambda_times, lambdas);
  internal->setMu(mu_times, mus);
  internal->setPhi(phi_times, phis);
  internal->setDelta(delta_times, deltas);
  internal->setEta(eta_times, etas);

  // Rcout << " number of eta matrices: " << etas.size() << std::endl;
  // Rcout << " number of eta rows: "     << etas[0].size() << std::endl;
  // Rcout << " number of eta cols: "     << etas[0][0].size() << std::endl;


};

// tree and data

void TensorPhyloExternal::setTree(const std::string &aNewickTree) {
  Rcout << "Setting tree." << std::endl;
  internal->setTree(aNewickTree);
}

void TensorPhyloExternal::setData() {
  stop("NOT IMPLEMENTED.");
}

////////////////////////
// numerical settings //
////////////////////////

void TensorPhyloExternal::setNumberOfThreads(size_t nThreads) {
  internal->setNumberOfThreads(nThreads);
}

void TensorPhyloExternal::setInitialDeltaT(double initDeltaT) {
  internal->setInitialDeltaT(initDeltaT);
}

void TensorPhyloExternal::setLikelihoodApproximator(int approxVersion) {
  internal->setLikelihoodApproximator( (approximatorVersion_t)approxVersion );
}

void TensorPhyloExternal::setSeed(size_t aSeed) {
  internal->setSeed(aSeed);
}

///////////////
// debugging //
///////////////

void TensorPhyloExternal::setDebugMode(int m) {
  internal->setDebugMode( (debugMode_t)m );
}

void TensorPhyloExternal::setSyncMonitors(const std::vector< double > &synchMonitoring) {
  internal->setSyncMonitors(synchMonitoring);
}

////////////////////
// model settings //
////////////////////

void TensorPhyloExternal::setApplyTreeLikCorrection(bool doApply) {
  internal->setApplyTreeLikCorrection(doApply);
}

void TensorPhyloExternal::setConditionalProbCompatibilityMode(bool setActive) {
  internal->setConditionalProbCompatibilityMode(setActive);
}


void TensorPhyloExternal::setConditionalProbabilityType(int condProb) {
  internal->setConditionalProbabilityType( (conditionalProbability_t)condProb );
}

////////////////
// likelihood //
////////////////

double TensorPhyloExternal::computeLogLikelihood() {

  // TODO: check for updates

  return internal->computeLogLikelihood();

}

////////////////
// converters //
////////////////

stdVectorXd TensorPhyloExternal::EigenToStd(VectorXd eig_vec) {

  // create the vector
  stdVectorXd res(eig_vec.data(), eig_vec.data() + eig_vec.size());

  return res;

}

VectorXd TensorPhyloExternal::StdToEigen(stdVectorXd std_vec) {

  // create the vector
  VectorXd res = Map<VectorXd, Eigen::Unaligned>(std_vec.data(), std_vec.size());

  return res;

}

stdMatrixXd TensorPhyloExternal::EigenToStd(MatrixXd eig_mat) {

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

MatrixXd TensorPhyloExternal::StdToEigen(stdMatrixXd std_mat) {

  // create the matrix
  MatrixXd mat( std_mat.size(), std_mat[0].size() );
  for(size_t r = 0; r < mat.rows(); ++r) {
    for(size_t c = 0; c < mat.cols(); ++c) {
      mat(r,c) = std_mat.at(r).at(c);
    }
  }

  return mat;

}


std::vector<stdMatrixXd> TensorPhyloExternal::ArmaToStd(arma::cube arma_cube) {

  // create the vector
  std::vector< stdMatrixXd > arr(arma_cube.n_slices);

  // loop over dimensions of the cube
  for(size_t i = 0; i < arma_cube.n_slices; ++i) {

    // get the slice (matrix)
    arma::mat rate_matrix = arma_cube.slice(i);

    // create the matrix
    stdMatrixXd std_mat;
    for(size_t r = 0; r < rate_matrix.n_rows; ++r) {
      stdVectorXd row;
      for(size_t c = 0; c < rate_matrix.n_cols; ++c) {
        row.push_back( rate_matrix(r, c) );
      }
      std_mat.push_back(row);
    }

    // push it back
    arr.at(i) = std_mat;

  }

  return arr;

}

arma::cube TensorPhyloExternal::StdToArma(std::vector<stdMatrixXd> std_cube) {

  // initialize the array
  arma::cube arr( std_cube[0].size(), std_cube[0][0].size(), std_cube.size());

  // fill in the values
  for(size_t t = 0; t < arr.n_slices; ++t) {
    for(size_t r = 0; r < arr.n_rows; ++r) {
      for(size_t c = 0; c < arr.n_cols; ++c) {
        arr(r, c, t) = std_cube[t][r][c];
      }
    }
  }

  // return the array
  return arr;

}

////////////////
// root prior //
////////////////

VectorXd TensorPhyloExternal::getRootPrior() {
  return StdToEigen(root_frequency);
}

void TensorPhyloExternal::setRootPrior(VectorXd new_root_freq) {

  // make sure root frequencies sum to 1
  if ( fabs( new_root_freq.lpNorm<1>() - 1.0) > 1e-10 ) {
    stop("Error setting root frequency. Frequencies must sum to 1.");
  }

  // check the number of states
  if ( new_root_freq.size() != dim ) {
    stop("Error setting root frequency. Number of frequencies should be equal to the number of states.");
  }

  // set the value
  internal->setRootPrior( EigenToStd(new_root_freq) );

}

////////////
// lambda //
////////////

MatrixXd TensorPhyloExternal::getLambda() {
  return StdToEigen(lambdas);
}

VectorXd TensorPhyloExternal::getLambdaTimes() {
  return StdToEigen(lambda_times);
}

void TensorPhyloExternal::setLambdaConstant(double new_lambda) {

  // set lambda times to empty
  lambda_times = stdVectorXd();

  // repeat the lambdas
  lambdas = stdMatrixXd(1, stdVectorXd(dim, new_lambda));

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set state varying
void TensorPhyloExternal::setLambdaStateVarying(VectorXd new_lambda) {

  // error handling
  // make sure the length of the vector matches the data
  if (new_lambda.size() != dim) {
    stop("Error setting speciation rates. Number of rates must equal the number of states.");
  }

  // set lambda times to empty
  lambda_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_lambda;

  lambdas = EigenToStd(ll);

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set time-varying lambda
void TensorPhyloExternal::setLambdaTimeVarying(VectorXd new_lambda_times, VectorXd new_lambda) {

  // error handling

  // make sure the number of times and states match
  if ( new_lambda.size() != (new_lambda_times.size() + 1) ) {
    stop("Error setting speciation rates. Number of change times must be 1 less than the number of speciation rates.");
  }

  // set the time variable
  lambda_times = EigenToStd(new_lambda_times);

  // copy the time-varying rates per state
  lambdas.resize( new_lambda.size());
  for(size_t i = 0; i < new_lambda.size(); ++i) {
    lambdas.at(i) = stdVectorXd(dim, new_lambda(i));
  }

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set time/state varying lambda
void TensorPhyloExternal::setLambdaTimeStateVarying(VectorXd new_lambda_times, MatrixXd new_lambda) {

  // do some error handling

  // make sure the number of times and states match
  if ( new_lambda.rows() != (new_lambda_times.size() + 1) ) {
    stop("Error setting speciation rates. Number of change times must be 1 less than the number of speciation rate vectors.");
  }

  // make sure the number of columns is correct
  if ( new_lambda.cols() != dim ) {
    stop("Error setting speciation rates. Number of rates per vector must equal the number of states.");
  }

  // set the time variable
  lambda_times = EigenToStd(new_lambda_times);

  // set the rate variable
  lambdas = EigenToStd(new_lambda);

  // set value
  internal->setLambda(lambda_times, lambdas);

}



////////
// mu //
////////

MatrixXd TensorPhyloExternal::getMu() {
  return StdToEigen(mus);
}

VectorXd TensorPhyloExternal::getMuTimes() {
  return StdToEigen(mu_times);
}

void TensorPhyloExternal::setMuConstant(double new_mu) {

  // set mu times to empty
  mu_times = stdVectorXd();

  // repeat the mus
  mus = stdMatrixXd(1, stdVectorXd(dim, new_mu));

  // set the value
  internal->setMu(mu_times, mus);

}

// set state varying
void TensorPhyloExternal::setMuStateVarying(VectorXd new_mu) {

  // error handling
  // make sure the length of the vector matches the data
  if (new_mu.size() != dim) {
    stop("Error setting extinction rates. Number of rates must equal the number of states.");
  }

  // set mu times to empty
  mu_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_mu;

  mus = EigenToStd(ll);

  // set the value
  internal->setMu(mu_times, mus);

}

// set time-varying mu
void TensorPhyloExternal::setMuTimeVarying(VectorXd new_mu_times, VectorXd new_mu) {

  // error handling

  // make sure the number of times and states match
  if ( new_mu.size() != (new_mu_times.size() + 1) ) {
    stop("Error setting extinction rates. Number of change times must be 1 less than the number of extinction rates.");
  }

  // set the time variable
  mu_times = EigenToStd(new_mu_times);

  // copy the time-varying rates per state
  mus.resize( new_mu.size());
  for(size_t i = 0; i < new_mu.size(); ++i) {
    mus.at(i) = stdVectorXd(dim, new_mu(i));
  }

  // set the value
  internal->setMu(mu_times, mus);

}

// set time/state varying mu
void TensorPhyloExternal::setMuTimeStateVarying(VectorXd new_mu_times, MatrixXd new_mu) {

  // do some error handling

  // make sure the number of times and states match
  if ( new_mu.rows() != (new_mu_times.size() + 1) ) {
    stop("Error setting extinction rates. Number of change times must be 1 less than the number of extinction rate vectors.");
  }

  // make sure the number of columns is correct
  if ( new_mu.cols() != dim ) {
    stop("Error setting extinction rates. Number of rates per vector must equal the number of states.");
  }

  // set the time variable
  mu_times = EigenToStd(new_mu_times);

  // set the rate variable
  mus = EigenToStd(new_mu);

  // set value
  internal->setMu(mu_times, mus);

}




/////////
// phi //
/////////

MatrixXd TensorPhyloExternal::getPhi() {
  return StdToEigen(phis);
}

VectorXd TensorPhyloExternal::getPhiTimes() {
  return StdToEigen(phi_times);
}

void TensorPhyloExternal::setPhiConstant(double new_phi) {

  // set phi times to empty
  phi_times = stdVectorXd();

  // repeat the phis
  phis = stdMatrixXd(1, stdVectorXd(dim, new_phi));

  // set the value
  internal->setPhi(phi_times, phis);

}

// set state varying
void TensorPhyloExternal::setPhiStateVarying(VectorXd new_phi) {

  // error handling
  // make sure the length of the vector matches the data
  if (new_phi.size() != dim) {
    stop("Error setting sampling rates. Number of rates must equal the number of states.");
  }

  // set phi times to empty
  phi_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_phi;

  phis = EigenToStd(ll);

  // set the value
  internal->setPhi(phi_times, phis);

}

// set time-varying phi
void TensorPhyloExternal::setPhiTimeVarying(VectorXd new_phi_times, VectorXd new_phi) {

  // error handling

  // make sure the number of times and states match
  if ( new_phi.size() != (new_phi_times.size() + 1) ) {
    stop("Error setting sampling rates. Number of change times must be 1 less than the number of sampling rates.");
  }

  // set the time variable
  phi_times = EigenToStd(new_phi_times);

  // copy the time-varying rates per state
  phis.resize( new_phi.size());
  for(size_t i = 0; i < new_phi.size(); ++i) {
    phis.at(i) = stdVectorXd(dim, new_phi(i));
  }

  // set the value
  internal->setPhi(phi_times, phis);

}

// set time/state varying phi
void TensorPhyloExternal::setPhiTimeStateVarying(VectorXd new_phi_times, MatrixXd new_phi) {

  // do some error handling

  // make sure the number of times and states match
  if ( new_phi.rows() != (new_phi_times.size() + 1) ) {
    stop("Error setting sampling rates. Number of change times must be 1 less than the number of sampling rate vectors.");
  }

  // make sure the number of columns is correct
  if ( new_phi.cols() != dim ) {
    stop("Error setting sampling rates. Number of rates per vector must equal the number of states.");
  }

  // set the time variable
  phi_times = EigenToStd(new_phi_times);

  // set the rate variable
  phis = EigenToStd(new_phi);

  // set value
  internal->setPhi(phi_times, phis);

}






///////////
// delta //
///////////


MatrixXd TensorPhyloExternal::getDelta() {
  return StdToEigen(deltas);
}

VectorXd TensorPhyloExternal::getDeltaTimes() {
  return StdToEigen(delta_times);
}

void TensorPhyloExternal::setDeltaConstant(double new_delta) {

  // set delta times to empty
  delta_times = stdVectorXd();

  // repeat the deltas
  deltas = stdMatrixXd(1, stdVectorXd(dim, new_delta));

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set state varying
void TensorPhyloExternal::setDeltaStateVarying(VectorXd new_delta) {

  // error handling
  // make sure the length of the vector matches the data
  if (new_delta.size() != dim) {
    stop("Error setting destructive-sampling rates. Number of rates must equal the number of states.");
  }

  // set delta times to empty
  delta_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_delta;

  deltas = EigenToStd(ll);

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set time-varying delta
void TensorPhyloExternal::setDeltaTimeVarying(VectorXd new_delta_times, VectorXd new_delta) {

  // error handling

  // make sure the number of times and states match
  if ( new_delta.size() != (new_delta_times.size() + 1) ) {
    stop("Error setting destructive-sampling rates. Number of change times must be 1 less than the number of destructive-sampling rates.");
  }

  // set the time variable
  delta_times = EigenToStd(new_delta_times);

  // copy the time-varying rates per state
  deltas.resize( new_delta.size());
  for(size_t i = 0; i < new_delta.size(); ++i) {
    deltas.at(i) = stdVectorXd(dim, new_delta(i));
  }

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set time/state varying delta
void TensorPhyloExternal::setDeltaTimeStateVarying(VectorXd new_delta_times, MatrixXd new_delta) {

  // do some error handling

  // make sure the number of times and states match
  if ( new_delta.rows() != (new_delta_times.size() + 1) ) {
    stop("Error setting destructive-sampling rates. Number of change times must be 1 less than the number of destructive-sampling rate vectors.");
  }

  // make sure the number of columns is correct
  if ( new_delta.cols() != dim ) {
    stop("Error setting destructive-sampling rates. Number of rates per vector deltast equal the number of states.");
  }

  // set the time variable
  delta_times = EigenToStd(new_delta_times);

  // set the rate variable
  deltas = EigenToStd(new_delta);

  // set value
  internal->setDelta(delta_times, deltas);

}





/////////
// eta //
/////////

arma::cube TensorPhyloExternal::getEta() {
  return StdToArma(etas);
}

VectorXd TensorPhyloExternal::getEtaTimes() {
  return StdToEigen(eta_times);
}

void TensorPhyloExternal::setEtaConstantEqual(double new_eta) {

  // set delta times to empty
  eta_times = stdVectorXd();

  // create a matrix
  MatrixXd tmp       = MatrixXd::Constant(dim, dim, new_eta);
  tmp.diagonal()     = ((double)dim - 1) * VectorXd::Constant(dim, -new_eta);
  stdMatrixXd tmpStd = EigenToStd(tmp);

  // create the vector of etas
  etas = std::vector<stdMatrixXd>(1, tmpStd);

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaConstantUnequal(MatrixXd new_eta) {

  // check the dimensionality
  if ( new_eta.cols() != dim || new_eta.rows() != dim ) {
    stop("Error setting transition rates. Rate matrix must be X by X, where X is the number of states.");
  }

  // set delta times to empty
  eta_times = stdVectorXd();

  // create the matrix
  stdMatrixXd tmpStd = EigenToStd(new_eta);

  // create the vector of etas
  etas = std::vector<stdMatrixXd>(1, tmpStd);

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaTimeVaryingEqual(VectorXd new_eta_times, VectorXd new_eta) {

  // error handling

  // make sure the number of times and states match
  if ( new_eta.size() != (new_eta_times.size() + 1) ) {
    stop("Error setting transition rates. Number of change times must be 1 less than the number of transition rates.");
  }

  // set the time variable
  eta_times = EigenToStd(new_eta_times);

  // make a rate matrix for each time
  etas.resize( new_eta.size() );
  for(size_t i = 0; i < new_eta.size(); ++i) {
    MatrixXd tmp   = MatrixXd::Constant(dim, dim, new_eta[i]);
    tmp.diagonal() = ((double)dim - 1) * VectorXd::Constant(dim, -new_eta[i]);
    etas.at(i) = EigenToStd(tmp);
  }

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaTimeVaryingUnequal(VectorXd new_eta_times, arma::cube new_eta) {

  // error handling

  // make sure the number of times and states match
  if ( new_eta.n_slices != (new_eta_times.size() + 1) ) {
    stop("Error setting transition rates. Number of change times must be 1 less than the number of transition rates.");
  }

  // check the dimensionality
  if ( new_eta.n_rows != dim || new_eta.n_cols != dim ) {
    stop("Error setting transition rates. Rate matrix must be X by X, where X is the number of states.");
  }

  // set the time variable
  eta_times = EigenToStd(new_eta_times);

  // create the "std" cube
  std::vector<stdMatrixXd> tmp = ArmaToStd(new_eta);
  etas = tmp;

  // set the value
  internal->setEta(eta_times, etas);

}





// ///////////
// // omega //
// ///////////
//
// void TensorPhyloExternal::setOmega(size_t aNState, const NumericVector &times, const std::vector< eventMap_t > &omegas) {
//   // internal->setOmega(aNState, times, omegas);
// }
//
// /////////////////////
// // mass speciation //
// /////////////////////
//
//
// void TensorPhyloExternal::setMassSpeciationEvents(const NumericVector &massSpeciationTimes, const NumericMatrix &massSpeciationProb) {
//   // internal->setMassSpeciationEvents(massSpeciationTimes, massSpeciationProb);
// }
//
// // mass extinction
//
// void TensorPhyloExternal::setMassExtinctionEvents(const NumericVector &massExtinctionTimes, const NumericMatrix &massExtinctionProb) {
//   // internal->setMassExtinctionEvents(massExtinctionTimes, massExtinctionProb);
// }
//
// // mass-extinction-induced state change
//
// void TensorPhyloExternal::setMassExtinctionStateChangeProb(const std::vector< NumericMatrix> &massExtinctionStateChangeProb) {
//   // internal->setMassExtinctionStateChangeProb(massExtinctionStateChangeProb);
// }
//
// // mass-sampling events
//
// void TensorPhyloExternal::setMassSamplingEvents(const NumericVector &massSamplingTimes, const NumericMatrix &massSamplingProb) {
//   // internal->setMassSamplingEvents(massSamplingTimes, massSamplingProb);
// }
//
// // mass-destructive-sampling events
//
// void TensorPhyloExternal::setMassDestrSamplingEvents(const NumericVector &massDestrSamplingTimes, const NumericMatrix &massDestrSamplingProb) {
//   // internal->setMassDestrSamplingEvents(massDestrSamplingTimes, massDestrSamplingProb);
// }
















// END: TensorPhyloExternal.cpp
