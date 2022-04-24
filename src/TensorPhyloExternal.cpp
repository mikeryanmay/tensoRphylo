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
  rf_dirty(true),
  rootFrequency( VectorXd::Ones(dim) / dim ),
  lambda_dirty(true),
  lambdas( MatrixXd::Ones(dim, 1) ),
  lambda_times(0),
  mu_dirty(true),
  mus( MatrixXd::Zero(dim, 1) ),
  mu_times(0),
  phi_dirty(true),
  phis( MatrixXd::Zero(dim, 1) ),
  phi_times(0),
  delta_dirty(true),
  deltas( MatrixXd::Zero(dim, 1) ),
  delta_times(0),
  eta_dirty(true),
  etas( arma::cube(dim, dim, 1, arma::fill::zeros) ),
  eta_times(0)
{
  // create the internal
  internal = DistributionHandlerImpl::create();
};

// tree and data

void TensorPhyloExternal::setTree(const std::string &aNewickTree) {
  Rcout << "Setting tree." << std::endl;
  internal->setTree(aNewickTree);
}

void TensorPhyloExternal::setData() {
  stop("NOT IMPLEMENTED.");
}

// numerical settings

void TensorPhyloExternal::setNumberOfThreads(size_t nThreads) {
  internal->setNumberOfThreads(nThreads);
}

void TensorPhyloExternal::setInitialDeltaT(double initDeltaT) {
  internal->setInitialDeltaT(initDeltaT);
}

// debugging

void TensorPhyloExternal::setDebugMode(int m) {
  internal->setDebugMode( (debugMode_t)m );
}

void TensorPhyloExternal::setSyncMonitors(const std::vector< double > &synchMonitoring) {
  internal->setSyncMonitors(synchMonitoring);
}

void TensorPhyloExternal::setSeed(size_t aSeed) {
  internal->setSeed(aSeed);
}

// model settings
void TensorPhyloExternal::setApplyTreeLikCorrection(bool doApply) {
  internal->setApplyTreeLikCorrection(doApply);
}

void TensorPhyloExternal::setConditionalProbCompatibilityMode(bool setActive) {
  internal->setConditionalProbCompatibilityMode(setActive);
}

void TensorPhyloExternal::setLikelihoodApproximator(int approxVersion) {
  internal->setLikelihoodApproximator( (approximatorVersion_t)approxVersion );
}

void TensorPhyloExternal::setConditionalProbabilityType(int condProb) {
  internal->setConditionalProbabilityType( (conditionalProbability_t)condProb );
}

////////////////
// converters //
////////////////

stdVectorXd TensorPhyloExternal::EigenToStd(VectorXd eig_vec) {

  // create the vector
  stdVectorXd res(eig_vec.data(), eig_vec.data() + eig_vec.size());

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
    arr.push_back(std_mat);

  }

  return arr;

}

////////////////
// root prior //
////////////////

VectorXd TensorPhyloExternal::getRootPrior() {
  return rootFrequency;
}

void TensorPhyloExternal::setRootPrior(VectorXd new_root_freq) {
  rf_dirty = true;
  rootFrequency = new_root_freq;
}

void TensorPhyloExternal::updateRootPrior() {

  // only update if something has changed
  if ( rf_dirty == false ) {
    return;
  }

  // check the number of states
  if ( rootFrequency.size() != dim ) {
    stop("Error setting root frequency. Number of frequencies should be equal to the number of states.");
  }

  // check that we sum to 1
  if ( fabs(rootFrequency.lpNorm<1>() - 1.0) > 1e-10 ) {
    stop("Error setting root frequency. Frequencies must sum to 1.");
  }

  // convert times to std::vector
  stdVectorXd pi = EigenToStd(rootFrequency);

  // set the values
  internal->setRootPrior(pi);

  // reset the dirty flag
  rf_dirty = false;

}

////////////
// lambda //
////////////

MatrixXd TensorPhyloExternal::getLambda() {
  return lambdas;
}

void TensorPhyloExternal::setLambda(MatrixXd new_lambda) {
  lambda_dirty = true;
  lambdas = new_lambda;
}

VectorXd TensorPhyloExternal::getLambdaTimes() {
  return lambda_times;
}

void TensorPhyloExternal::setLambdaTimes(VectorXd new_lambda_times) {
  lambda_dirty = true;
  lambda_times  = new_lambda_times;
}

void TensorPhyloExternal::updateLambdas() {

  // only update if something has changed
  if ( lambda_dirty == false ) {
    return;
  }

  // check the number of times
  size_t num_times = lambda_times.size();
  if ( lambdas.cols() != (num_times + 1) ) {
    stop("Error setting speciation rates. Number of change times must be 1 less than the number of columns in the speciation matrix.");
  }

  // check the number of states
  if ( lambdas.rows() != dim ) {
    stop("Error setting speciation rates. Number of rates should (rows) should be equal to the number of states.");
  }

  // convert times to std::vector
  stdVectorXd lt = EigenToStd(lambda_times);

  // convert rates to matrix
  stdMatrixXd ls = EigenToStd(lambdas);

  // set the values
  internal->setLambda(lt, ls);

  // reset the dirty flag
  lambda_dirty = false;

}

////////
// mu //
////////

MatrixXd TensorPhyloExternal::getMu() {
  return mus;
}

void TensorPhyloExternal::setMu(MatrixXd new_mu) {
  mu_dirty = true;
  mus = new_mu;
}

VectorXd TensorPhyloExternal::getMuTimes() {
  return mu_times;
}

void TensorPhyloExternal::setMuTimes(VectorXd new_mu_times) {
  mu_dirty = true;
  mu_times  = new_mu_times;
}

void TensorPhyloExternal::updateMus() {

  // only update if something has changed
  if ( mu_dirty == false ) {
    return;
  }

  // check the number of times
  size_t num_times = mu_times.size();
  if ( mus.cols() != (num_times + 1) ) {
    stop("Error setting extinction rates. Number of change times must be 1 less than the number of columns in the speciation matrix.");
  }

  // check the number of states
  if ( mus.rows() != dim ) {
    stop("Error setting extinction rates. Number of rates should (rows) should be equal to the number of states.");
  }

  // convert times to std::vector
  stdVectorXd mt = EigenToStd(mu_times);

  // convert rates to matrix
  stdMatrixXd ms = EigenToStd(mus);

  // set the values
  internal->setMu(mt, ms);

  // reset the dirty flag
  mu_dirty = false;

}


/////////
// phi //
/////////

MatrixXd TensorPhyloExternal::getPhi() {
  return phis;
}

void TensorPhyloExternal::setPhi(MatrixXd new_phi) {
  phi_dirty = true;
  phis = new_phi;
}

VectorXd TensorPhyloExternal::getPhiTimes() {
  return phi_times;
}

void TensorPhyloExternal::setPhiTimes(VectorXd new_phi_times) {
  phi_dirty = true;
  phi_times  = new_phi_times;
}

void TensorPhyloExternal::updatePhis() {

  // only update if something has changed
  if ( phi_dirty == false ) {
    return;
  }

  // check the number of times
  size_t num_times = phi_times.size();
  if ( phis.cols() != (num_times + 1) ) {
    stop("Error setting sampling rates. Number of change times must be 1 less than the number of columns in the speciation matrix.");
  }

  // check the number of states
  if ( phis.rows() != dim ) {
    stop("Error setting sampling rates. Number of rates should (rows) should be equal to the number of states.");
  }

  // convert times to std::vector
  stdVectorXd pt = EigenToStd(phi_times);

  // convert rates to matrix
  stdMatrixXd ps = EigenToStd(phis);

  // set the values
  internal->setPhi(pt, ps);

  // reset the dirty flag
  phi_dirty = false;

}


///////////
// delta //
///////////

MatrixXd TensorPhyloExternal::getDelta() {
  return deltas;
}

void TensorPhyloExternal::setDelta(MatrixXd new_delta) {
  delta_dirty = true;
  deltas = new_delta;
}

VectorXd TensorPhyloExternal::getDeltaTimes() {
  return delta_times;
}

void TensorPhyloExternal::setDeltaTimes(VectorXd new_delta_times) {
  delta_dirty = true;
  delta_times  = new_delta_times;
}

void TensorPhyloExternal::updateDeltas() {

  // only update if something has changed
  if ( delta_dirty == false ) {
    return;
  }

  // check the number of times
  size_t num_times = delta_times.size();
  if ( deltas.cols() != (num_times + 1) ) {
    stop("Error setting destructive-sampling rates. Number of change times must be 1 less than the number of columns in the speciation matrix.");
  }

  // check the number of states
  if ( deltas.rows() != dim ) {
    stop("Error setting destructive-sampling rates. Number of rates should (rows) should be equal to the number of states.");
  }

  // convert times to std::vector
  stdVectorXd dt = EigenToStd(delta_times);

  // convert rates to matrix
  stdMatrixXd ds = EigenToStd(deltas);

  // set the values
  internal->setDelta(dt, ds);

  // reset the dirty flag
  delta_dirty = false;

}

/////////
// eta //
/////////

arma::cube TensorPhyloExternal::getEta() {
  return etas;
}

void TensorPhyloExternal::setEta(arma::cube new_eta) {
  eta_dirty = true;
  etas = new_eta;
}

VectorXd TensorPhyloExternal::getEtaTimes() {
  return eta_times;
}

void TensorPhyloExternal::setEtaTimes(VectorXd new_eta_times) {
  eta_dirty = true;
  eta_times  = new_eta_times;
}

void TensorPhyloExternal::updateEtas() {

  // only update if something has changed
  if ( eta_dirty == false ) {
    return;
  }

  // check the number of times
  size_t num_times = eta_times.size();
  if ( etas.n_slices != (num_times + 1) ) {
    stop("Error setting transition rates. Number of rate matrices must be 1 less than the number of columns in the speciation matrix.");
  }

  // // check the number of states
  if ( etas.n_rows != dim || etas.n_cols != dim ) {
    stop("Error setting transition rates. Number of rates should (rows and columbns) should be equal to the number of states.");
  }

  // convert times to std::vector
  stdVectorXd et = EigenToStd(eta_times);

  // convert transition rates to std::vector< stdMatrixXd >
  std::vector<stdMatrixXd> qv = ArmaToStd(etas);

  // // set the values
  internal->setEta(et, qv);

  // reset the dirty flag
  eta_dirty = false;

}


///////////
// omega //
///////////

void TensorPhyloExternal::setOmega(size_t aNState, const NumericVector &times, const std::vector< eventMap_t > &omegas) {
  // internal->setOmega(aNState, times, omegas);
}

/////////////////////
// mass speciation //
/////////////////////


void TensorPhyloExternal::setMassSpeciationEvents(const NumericVector &massSpeciationTimes, const NumericMatrix &massSpeciationProb) {
  // internal->setMassSpeciationEvents(massSpeciationTimes, massSpeciationProb);
}

// mass extinction

void TensorPhyloExternal::setMassExtinctionEvents(const NumericVector &massExtinctionTimes, const NumericMatrix &massExtinctionProb) {
  // internal->setMassExtinctionEvents(massExtinctionTimes, massExtinctionProb);
}

// mass-extinction-induced state change

void TensorPhyloExternal::setMassExtinctionStateChangeProb(const std::vector< NumericMatrix> &massExtinctionStateChangeProb) {
  // internal->setMassExtinctionStateChangeProb(massExtinctionStateChangeProb);
}

// mass-sampling events

void TensorPhyloExternal::setMassSamplingEvents(const NumericVector &massSamplingTimes, const NumericMatrix &massSamplingProb) {
  // internal->setMassSamplingEvents(massSamplingTimes, massSamplingProb);
}

// mass-destructive-sampling events

void TensorPhyloExternal::setMassDestrSamplingEvents(const NumericVector &massDestrSamplingTimes, const NumericMatrix &massDestrSamplingProb) {
  // internal->setMassDestrSamplingEvents(massDestrSamplingTimes, massDestrSamplingProb);
}

// likelihoods
double TensorPhyloExternal::computeLogLikelihood() {

  // TODO: check for updates

  return internal->computeLogLikelihood();
}















// END: TensorPhyloExternal.cpp
