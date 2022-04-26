#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "TensorPhyloExternal.h"
#include "TensorPhyloUtils.h"

using namespace Rcpp;
// using namespace Eigen;
using namespace TensorPhylo::Interface;

// constructors/destructors //
TensorPhyloExternal::TensorPhyloExternal(size_t dim_) :
  safe(true),
  dim(dim_)
{
  // create the internal
  internal = DistributionHandlerImpl::create();

  // TODO: use setters to define defaults
  // so we can use validators

  // default root frequency is flat
  setRootPriorFlat();

  // default speciation rate is 1.0
  setLambdaConstant(1.0);

  // default extinction rate is 0.0
  setMuConstant(1.0);

  // default sampling rates are zero
  setPhiConstant(0.0);
  setDeltaConstant(0.0);

  // default transition rate is 1.0
  setEtaConstantEqual(1.0);

  // default mass-sampling probability is 1 at the present
  setRhoPresent(1.0);

  // parameters not set (not in default model):
  // omega
  // upsilon
  // gamma nor zeta
  // xi

};

///////////////////
// tree and data //
///////////////////

// TODO: this stuff should be private, and the
// tree and data should be provided to the ctor

void TensorPhyloExternal::setTree(const std::string &aNewickTree) {
  Rcout << "Setting tree." << std::endl;
  internal->setTree(aNewickTree);
}

void TensorPhyloExternal::setData() {
  stop("NOT YET IMPLEMENTED.");
}

////////////////////////
// numerical settings //
////////////////////////

void TensorPhyloExternal::setSafeMode(bool safe_mode) {
  safe = safe_mode;
}

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
// root prior //
////////////////

VectorXd TensorPhyloExternal::getRootPrior() {
  return TensorPhyloUtils::StdToEigen(root_frequency);
}

void TensorPhyloExternal::setRootPriorFlat() {

  // make a flat vector
  root_frequency = stdVectorXd(dim, 1.0 / (double)dim);

  // update
  internal->setRootPrior(root_frequency);

}

void TensorPhyloExternal::setRootPrior(VectorXd new_root_freq) {

  if ( safe ) {

    // check that this is a valid simplex
    if ( TensorPhyloUtils::isSimplex(new_root_freq) == false ) {
      stop("Error setting root frequency. Frequencies must sum to 1.");
    }

    // check the number of states
    if ( TensorPhyloUtils::hasDimensions(new_root_freq, dim) == false ) {
      stop("Error setting root frequency. Number of frequencies should be equal to the number of states.");
    }

  }

  root_frequency = TensorPhyloUtils::EigenToStd(new_root_freq);

  // set the value
  internal->setRootPrior(root_frequency);

}

////////////
// lambda //
////////////

MatrixXd TensorPhyloExternal::getLambda() {
  return TensorPhyloUtils::StdToEigen(lambdas);
}

VectorXd TensorPhyloExternal::getLambdaTimes() {
  return TensorPhyloUtils::StdToEigen(lambda_times);
}

void TensorPhyloExternal::setLambdaConstant(double new_lambda) {

  if ( safe ) {

    // check that time is valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda) == false ) {
      stop("Error setting speciation rates. Rates must be strictly non-negative.");
    }

  }

  // set lambda times to empty
  lambda_times = stdVectorXd();

  // repeat the lambdas
  lambdas = stdMatrixXd(1, stdVectorXd(dim, new_lambda));

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set state varying
void TensorPhyloExternal::setLambdaStateVarying(VectorXd new_lambda) {

  if ( safe ) {

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda) == false ) {
      stop("Error setting speciation rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_lambda, dim) == false ) {
      stop("Error setting speciation rates. Number of rates must equal the number of states.");
    }

  }

  // set lambda times to empty
  lambda_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_lambda;

  lambdas = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set time-varying lambda
void TensorPhyloExternal::setLambdaTimeVarying(VectorXd new_lambda_times, VectorXd new_lambda) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda_times) == false ) {
      stop("Error setting speciation rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda) == false ) {
      stop("Error setting speciation rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_lambda, new_lambda_times.size() + 1) == false ) {
      stop("Error setting speciation rates. Number of change times must be 1 less than the number of speciation rates.");
    }

  }

  // set the time variable
  lambda_times = TensorPhyloUtils::EigenToStd(new_lambda_times);

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

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda_times) == false ) {
      stop("Error setting speciation rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda) == false ) {
      stop("Error setting speciation rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_lambda, new_lambda_times.size() + 1, dim) == false ) {
      stop("Error setting speciation rates. Number of rows must be equal to the number of change times plus 1, and the number of columns must be equal to the number of states.");
    }

  }

  // set the time variable
  lambda_times = TensorPhyloUtils::EigenToStd(new_lambda_times);

  // set the rate variable
  lambdas = TensorPhyloUtils::EigenToStd(new_lambda);

  // set value
  internal->setLambda(lambda_times, lambdas);

}



////////
// mu //
////////

MatrixXd TensorPhyloExternal::getMu() {
  return TensorPhyloUtils::StdToEigen(mus);
}

VectorXd TensorPhyloExternal::getMuTimes() {
  return TensorPhyloUtils::StdToEigen(mu_times);
}

void TensorPhyloExternal::setMuConstant(double new_mu) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu) == false ) {
      stop("Error setting extinction rates. Rates must be strictly non-negative.");
    }

  }

  // set mu times to empty
  mu_times = stdVectorXd();

  // repeat the mus
  mus = stdMatrixXd(1, stdVectorXd(dim, new_mu));

  // set the value
  internal->setMu(mu_times, mus);

}

// set state varying
void TensorPhyloExternal::setMuStateVarying(VectorXd new_mu) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu) == false ) {
      stop("Error setting extinction rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_mu, dim) == false ) {
      stop("Error setting extinction rates. Number of rates must equal the number of states.");
    }

  }

  // set mu times to empty
  mu_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_mu;

  mus = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setMu(mu_times, mus);

}

// set time-varying mu
void TensorPhyloExternal::setMuTimeVarying(VectorXd new_mu_times, VectorXd new_mu) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu_times) == false ) {
      stop("Error setting extinction rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu) == false ) {
      stop("Error setting extinction rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_mu, new_mu_times.size() + 1) == false ) {
      stop("Error setting extinction rates. Number of change times must be 1 less than the number of extinction rates.");
    }

  }

  // set the time variable
  mu_times = TensorPhyloUtils::EigenToStd(new_mu_times);

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

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu_times) == false ) {
      stop("Error setting speciation rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu) == false ) {
      stop("Error setting extinction rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_mu, new_mu_times.size() + 1, dim) == false ) {
      stop("Error setting extinction rates. Number of rows must be equal to the number of change times plus 1, and the number of columns must be equal to the number of states.");
    }

  }

  // set the time variable
  mu_times = TensorPhyloUtils::EigenToStd(new_mu_times);

  // set the rate variable
  mus = TensorPhyloUtils::EigenToStd(new_mu);

  // set value
  internal->setMu(mu_times, mus);

}




/////////
// phi //
/////////

MatrixXd TensorPhyloExternal::getPhi() {
  return TensorPhyloUtils::StdToEigen(phis);
}

VectorXd TensorPhyloExternal::getPhiTimes() {
  return TensorPhyloUtils::StdToEigen(phi_times);
}

void TensorPhyloExternal::setPhiConstant(double new_phi) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
      stop("Error setting sampling rates. Rates must be strictly non-negative.");
    }

  }

  // set lambda times to empty
  phi_times = stdVectorXd();

  // repeat the phis
  phis = stdMatrixXd(1, stdVectorXd(dim, new_phi));

  // set the value
  internal->setPhi(phi_times, phis);

}

// set state varying
void TensorPhyloExternal::setPhiStateVarying(VectorXd new_phi) {

  if ( safe ) {

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
      stop("Error setting sampling rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_phi, dim) == false ) {
      Rcout << new_phi.size() << " -- " << dim << std::endl;
      stop("Error setting sampling rates. Number of rates must equal the number of states.");
    }

  }

  // set phi times to empty
  phi_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_phi;

  phis = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setPhi(phi_times, phis);

}

// set time-varying phi
void TensorPhyloExternal::setPhiTimeVarying(VectorXd new_phi_times, VectorXd new_phi) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi_times) == false ) {
      stop("Error setting sampling rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
      stop("Error setting sampling rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_phi, new_phi_times.size() + 1) == false ) {
      stop("Error setting sampling rates. Number of change times must be 1 less than the number of extinction rates.");
    }

  }

  // set the time variable
  phi_times = TensorPhyloUtils::EigenToStd(new_phi_times);

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

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi_times) == false ) {
      stop("Error setting sampling rates. Times must be strictly non-negative.");
    }

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
      stop("Error setting sampling rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_phi, new_phi_times.size() + 1, dim) == false ) {
      stop("Error setting sampling rates. Number of rows must be equal to the number of change times plus 1, and the number of columns must be equal to the number of states.");
    }

  }

  // set the time variable
  phi_times = TensorPhyloUtils::EigenToStd(new_phi_times);

  // set the rate variable
  phis = TensorPhyloUtils::EigenToStd(new_phi);

  // set value
  internal->setPhi(phi_times, phis);

}






///////////
// delta //
///////////


MatrixXd TensorPhyloExternal::getDelta() {
  return TensorPhyloUtils::StdToEigen(deltas);
}

VectorXd TensorPhyloExternal::getDeltaTimes() {
  return TensorPhyloUtils::StdToEigen(delta_times);
}

void TensorPhyloExternal::setDeltaConstant(double new_delta) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
      stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
    }

  }

  // set delta times to empty
  delta_times = stdVectorXd();

  // repeat the deltas
  deltas = stdMatrixXd(1, stdVectorXd(dim, new_delta));

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set state varying
void TensorPhyloExternal::setDeltaStateVarying(VectorXd new_delta) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
      stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_delta, dim) == false ) {
      Rcout << new_delta.size() << " -- " << dim << std::endl;
      stop("Error setting destructive-sampling rates. Number of rates must equal the number of states.");
    }

  }

  // set delta times to empty
  delta_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_delta;

  deltas = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set time-varying delta
void TensorPhyloExternal::setDeltaTimeVarying(VectorXd new_delta_times, VectorXd new_delta) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta_times) == false ) {
      stop("Error setting destructive-sampling rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
      stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_delta, new_delta_times.size() + 1) == false ) {
      stop("Error setting destructive-sampling rates. Number of change times must be 1 less than the number of extinction rates.");
    }

  }

  // set the time variable
  delta_times = TensorPhyloUtils::EigenToStd(new_delta_times);

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

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta_times) == false ) {
      stop("Error setting destructive-sampling rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
      stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_delta_times, new_delta.size() + 1, dim) == false ) {
      stop("Error setting destructive-sampling rates. Number of rows must be equal to the number of change times plus 1, and the number of columns must be equal to the number of states.");
    }

  }

  // set the time variable
  delta_times = TensorPhyloUtils::EigenToStd(new_delta_times);

  // set the rate variable
  deltas = TensorPhyloUtils::EigenToStd(new_delta);

  // set value
  internal->setDelta(delta_times, deltas);

}





/////////
// eta //
/////////

arma::cube TensorPhyloExternal::getEta() {
  return TensorPhyloUtils::StdToArma(etas);
}

VectorXd TensorPhyloExternal::getEtaTimes() {
  return TensorPhyloUtils::StdToEigen(eta_times);
}

void TensorPhyloExternal::setEtaConstantEqual(double new_eta) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta) == false ) {
      stop("Error setting transition rates. Rates must be strictly non-negative.");
    }

  }

  // set delta times to empty
  eta_times = stdVectorXd();

  // create a matrix
  MatrixXd tmp       = MatrixXd::Constant(dim, dim, new_eta);
  tmp.diagonal()     = ((double)dim - 1) * VectorXd::Constant(dim, -new_eta);
  stdMatrixXd tmpStd = TensorPhyloUtils::EigenToStd(tmp);

  // create the vector of etas
  etas = std::vector<stdMatrixXd>(1, tmpStd);

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaConstantUnequal(MatrixXd new_eta) {

  if ( safe ) {

    // make sure this is a rate matrix
    if ( TensorPhyloUtils::isRateMatrix(new_eta) == false ) {
      stop("Error setting transition rates. Rate matrix must be square, and the rows must sum to zero.");
    }

  }

  // set delta times to empty
  eta_times = stdVectorXd();

  // create the matrix
  stdMatrixXd tmpStd = TensorPhyloUtils::EigenToStd(new_eta);

  // create the vector of etas
  etas = std::vector<stdMatrixXd>(1, tmpStd);

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaTimeVaryingEqual(VectorXd new_eta_times, VectorXd new_eta) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta_times) == false ) {
      stop("Error setting transition rates. Times must be strictly non-negative.");
    }

    // check that rates are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta) == false ) {
      stop("Error setting transition rates. Rates must be strictly non-negative.");
    }

    // check the size
    if ( TensorPhyloUtils::hasDimensions(new_eta, new_eta_times.size() + 1) == false ) {
      stop("Error setting transition rates. Number of change times must be 1 less than the number of transition rates.");
    }

  }

  // set the time variable
  eta_times = TensorPhyloUtils::EigenToStd(new_eta_times);

  // make a rate matrix for each time
  etas.resize( new_eta.size() );
  for(size_t i = 0; i < new_eta.size(); ++i) {
    MatrixXd tmp   = MatrixXd::Constant(dim, dim, new_eta[i]);
    tmp.diagonal() = ((double)dim - 1) * VectorXd::Constant(dim, -new_eta[i]);
    etas.at(i) = TensorPhyloUtils::EigenToStd(tmp);
  }

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaTimeVaryingUnequal(VectorXd new_eta_times, arma::cube new_eta) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta_times) == false ) {
      stop("Error setting transition rates. Times must be strictly non-negative.");
    }

    // check dimensions of matrix
    if ( TensorPhyloUtils::hasDimensions(new_eta, dim, dim, new_eta_times.size() + 1) == false ) {
      stop("Error setting transition rates. Number of change times must be 1 less than the number of transition rates; each rate matrix must be square.");
    }

    // check each slice
    for(size_t i = 0; i < new_eta.n_slices; ++i) {
      if ( TensorPhyloUtils::isRateMatrix( new_eta.slice(i) ) == false ) {
        stop("Error setting transition rates. Rate matrix must be square, and the rows must sum to zero.");
      }
    }

  }

  // set the time variable
  eta_times = TensorPhyloUtils::EigenToStd(new_eta_times);

  // create the "std" cube
  etas = TensorPhyloUtils::ArmaToStd(new_eta);

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

/////////////////////
// mass speciation //
/////////////////////

MatrixXd TensorPhyloExternal::getUpsilon() {
  return TensorPhyloUtils::StdToEigen(upsilons);
}

VectorXd TensorPhyloExternal::getUpsilonTimes() {
  return TensorPhyloUtils::StdToEigen(upsilon_times);
}

void TensorPhyloExternal::setUpsilonConstant(VectorXd new_upsilon_times, VectorXd new_upsilon) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_upsilon_times) == false ) {
      stop("Error setting mass-speciation events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_upsilon, new_upsilon_times.size()) == false ) {
      stop("Error setting mass-speciation events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_upsilon) == false ) {
      stop("Error setting mass-speciation events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  upsilon_times = TensorPhyloUtils::EigenToStd(new_upsilon_times);

  // copy the time-varying rates per state
  upsilons.resize( new_upsilon.size());
  for(size_t i = 0; i < new_upsilon.size(); ++i) {
    upsilons.at(i) = stdVectorXd(dim, new_upsilon(i));
  }

  // set value
  internal->setMassSpeciationEvents(upsilon_times, upsilons);

}

void TensorPhyloExternal::setUpsilonStateVarying(VectorXd new_upsilon_times, MatrixXd new_upsilon) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_upsilon_times) == false ) {
      stop("Error setting mass-speciation events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_upsilon, new_upsilon_times.size(), dim) == false ) {
      stop("Error setting mass-speciation events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_upsilon) == false ) {
      stop("Error setting mass-speciation events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  upsilon_times = TensorPhyloUtils::EigenToStd(new_upsilon_times);

  // set the rate variable
  upsilons = TensorPhyloUtils::EigenToStd(new_upsilon);

  // set value
  internal->setMassSpeciationEvents(upsilon_times, upsilons);

}

/////////////////////
// mass extinction //
/////////////////////

MatrixXd TensorPhyloExternal::getGamma() {
  return TensorPhyloUtils::StdToEigen(gammas);
}

VectorXd TensorPhyloExternal::getGammaTimes() {
  return TensorPhyloUtils::StdToEigen(gamma_times);
}

void TensorPhyloExternal::setGammaConstant(VectorXd new_gamma_times, VectorXd new_gamma) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
      stop("Error setting mass-extinction events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_gamma, new_gamma_times.size()) == false ) {
      stop("Error setting mass-extinction events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_gamma) == false ) {
      stop("Error setting mass-extinction events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

  // copy the time-varying rates per state
  gammas.resize( new_gamma.size());
  for(size_t i = 0; i < new_gamma.size(); ++i) {
    gammas.at(i) = stdVectorXd(dim, new_gamma(i));
  }

  // set value
  internal->setMassSpeciationEvents(gamma_times, gammas);

}

void TensorPhyloExternal::setGammaStateVarying(VectorXd new_gamma_times, MatrixXd new_gamma) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
      stop("Error setting mass-extinction events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_gamma, new_gamma_times.size(), dim) == false ) {
      stop("Error setting mass-speciation events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_gamma) == false ) {
      stop("Error setting mass-speciation events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

  // set the rate variable
  gammas = TensorPhyloUtils::EigenToStd(new_gamma);

  // set value
  internal->setMassSpeciationEvents(gamma_times, gammas);

}

// TODO: mass-extinction-induced state change



///////////////////
// mass sampling //
///////////////////

MatrixXd TensorPhyloExternal::getRho() {
  return TensorPhyloUtils::StdToEigen(rhos);
}

VectorXd TensorPhyloExternal::getRhoTimes() {
  return TensorPhyloUtils::StdToEigen(rho_times);
}

void TensorPhyloExternal::setRhoPresent(double new_rho) {

  if ( safe ) {

    // check that rho is a probability
    if ( TensorPhyloUtils::isProbability(new_rho) == false ) {
      stop("Error setting mass-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set time to the present
  rho_times = stdVectorXd(1, 0);

  // create a matrix of sampling magnitudes
  MatrixXd tmp = MatrixXd::Constant(1, dim, new_rho);

  // set the matrix
  rhos = TensorPhyloUtils::EigenToStd(tmp);

  // set the value
  internal->setMassSamplingEvents(rho_times, rhos);

}

void TensorPhyloExternal::setRhoPresentStateVarying(VectorXd new_rho) {

  if ( safe ) {

    // check dimensions
    if ( TensorPhyloUtils::hasDimensions(new_rho, dim) == false ) {
      stop("Error setting mass-extinction events. Number of times must equal to the number of magnitudes.");
    }

    // check that rho is a probability
    if ( TensorPhyloUtils::isProbability(new_rho) == false ) {
      stop("Error setting mass-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set time to the present
  rho_times = stdVectorXd(1, 0);

  // create the matrix
  stdVectorXd rho_vec = TensorPhyloUtils::EigenToStd(new_rho);
  rhos = stdMatrixXd(1, rho_vec);

  // set the value
  internal->setMassSamplingEvents(rho_times, rhos);

}

void TensorPhyloExternal::setRhoConstant(VectorXd new_rho_times, VectorXd new_rho) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_rho_times) == false ) {
      stop("Error setting mass-sampling events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_rho, new_rho_times.size()) == false ) {
      stop("Error setting mass-sampling events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_rho) == false ) {
      stop("Error setting mass-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  rho_times = TensorPhyloUtils::EigenToStd(new_rho_times);

  // copy the time-varying rates per state
  rhos.resize( new_rho.size());
  for(size_t i = 0; i < new_rho.size(); ++i) {
    rhos.at(i) = stdVectorXd(dim, new_rho(i));
  }

  // set value
  internal->setMassSamplingEvents(rho_times, rhos);

}

void TensorPhyloExternal::setRhoStateVarying(VectorXd new_rho_times, MatrixXd new_rho) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_rho_times) == false ) {
      stop("Error setting mass-sampling events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_rho, new_rho_times.size(), dim) == false ) {
      stop("Error setting mass-sampling events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_rho) == false ) {
      stop("Error setting mass-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  rho_times = TensorPhyloUtils::EigenToStd(new_rho_times);

  // set the rate variable
  rhos = TensorPhyloUtils::EigenToStd(new_rho);

  // set value
  internal->setMassSamplingEvents(rho_times, rhos);

}


///////////////////////////////
// mass destructive-sampling //
///////////////////////////////

MatrixXd TensorPhyloExternal::getXi() {
  return TensorPhyloUtils::StdToEigen(xis);
}

VectorXd TensorPhyloExternal::getXiTimes() {
  return TensorPhyloUtils::StdToEigen(xi_times);
}

void TensorPhyloExternal::setXiConstant(VectorXd new_xi_times, VectorXd new_xi) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_xi_times) == false ) {
      stop("Error setting mass-destructive-sampling events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_xi, new_xi_times.size()) == false ) {
      stop("Error setting mass-destructive-sampling events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_xi) == false ) {
      stop("Error setting mass-destructive-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  xi_times = TensorPhyloUtils::EigenToStd(new_xi_times);

  // copy the time-varying rates per state
  xis.resize( new_xi.size());
  for(size_t i = 0; i < new_xi.size(); ++i) {
    xis.at(i) = stdVectorXd(dim, new_xi(i));
  }

  // set value
  internal->setMassDestrSamplingEvents(xi_times, xis);

}

void TensorPhyloExternal::setXiStateVarying(VectorXd new_xi_times, MatrixXd new_xi) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_xi_times) == false ) {
      stop("Error setting mass-destuctive-sampling events. Times must be strictly non-negative.");
    }

    // make sure the number of times and magnitudes match
    if ( TensorPhyloUtils::hasDimensions(new_xi, new_xi_times.size(), dim) == false ) {
      stop("Error setting mass-destuctive-sampling events. Number of times must equal to the number of magnitudes.");
    }

    // make sure the upsilons are probabilities
    if ( TensorPhyloUtils::isProbability(new_xi) == false ) {
      stop("Error setting mass-destuctive-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set the time variable
  xi_times = TensorPhyloUtils::EigenToStd(new_xi_times);

  // set the rate variable
  xis = TensorPhyloUtils::EigenToStd(new_xi);

  // set value
  internal->setMassDestrSamplingEvents(xi_times, xis);

}







// END: TensorPhyloExternal.cpp
