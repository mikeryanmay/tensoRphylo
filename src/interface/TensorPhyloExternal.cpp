#include <Rcpp.h>
#include <RcppEigen.h>
#include <typeinfo>
#include "TensorPhyloExternal.h"
#include "TensorPhyloUtils.h"
#include "CladoEvents.h"

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
  return internal->computeLogLikelihood();
}

////////////////
// root prior //
////////////////

void TensorPhyloExternal::setRootPriorFlat() {

  // make a flat vector
  stdVectorXd root_frequency = stdVectorXd(dim, 1.0 / (double)dim);

  // update
  internal->setRootPrior(root_frequency);

}

void TensorPhyloExternal::setRootPrior(const VectorXd& new_root_freq) {

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

  stdVectorXd root_frequency = TensorPhyloUtils::EigenToStd(new_root_freq);

  // set the value
  internal->setRootPrior(root_frequency);

}

////////////
// lambda //
////////////

void TensorPhyloExternal::setLambdaConstant(const double& new_lambda) {

  if ( safe ) {

    // check that time is valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda) == false ) {
      stop("Error setting speciation rates. Rates must be strictly non-negative.");
    }

  }

  // set lambda times to empty
  stdVectorXd lambda_times = stdVectorXd();

  // repeat the lambdas
  stdMatrixXd lambdas = stdMatrixXd(1, stdVectorXd(dim, new_lambda));

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set state varying
void TensorPhyloExternal::setLambdaStateVarying(const VectorXd& new_lambda) {

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
  stdVectorXd lambda_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_lambda;

  stdMatrixXd lambdas = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set time-varying lambda
void TensorPhyloExternal::setLambdaTimeVarying(const VectorXd& new_lambda_times, const VectorXd& new_lambda) {

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
  stdVectorXd lambda_times = TensorPhyloUtils::EigenToStd(new_lambda_times);

  // copy the time-varying rates per state
  stdMatrixXd lambdas(new_lambda.size());
  for(size_t i = 0; i < new_lambda.size(); ++i) {
    lambdas.at(i) = stdVectorXd(dim, new_lambda(i));
  }

  // set the value
  internal->setLambda(lambda_times, lambdas);

}

// set time/state varying lambda
void TensorPhyloExternal::setLambdaTimeStateVarying(const VectorXd& new_lambda_times, const MatrixXd& new_lambda) {

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
  stdVectorXd lambda_times = TensorPhyloUtils::EigenToStd(new_lambda_times);

  // set the rate variable
  stdMatrixXd lambdas = TensorPhyloUtils::EigenToStd(new_lambda);

  // set value
  internal->setLambda(lambda_times, lambdas);

}



////////
// mu //
////////

void TensorPhyloExternal::setMuConstant(const double& new_mu) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu) == false ) {
      stop("Error setting extinction rates. Rates must be strictly non-negative.");
    }

  }

  // set mu times to empty
  stdVectorXd mu_times = stdVectorXd();

  // repeat the mus
  stdMatrixXd mus = stdMatrixXd(1, stdVectorXd(dim, new_mu));

  // set the value
  internal->setMu(mu_times, mus);

}

// set state varying
void TensorPhyloExternal::setMuStateVarying(const VectorXd& new_mu) {

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
  stdVectorXd mu_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_mu;

  stdMatrixXd mus = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setMu(mu_times, mus);

}

// set time-varying mu
void TensorPhyloExternal::setMuTimeVarying(const VectorXd& new_mu_times, const VectorXd& new_mu) {

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
  stdVectorXd mu_times = TensorPhyloUtils::EigenToStd(new_mu_times);

  // copy the time-varying rates per state
  stdMatrixXd mus(new_mu.size());
  for(size_t i = 0; i < new_mu.size(); ++i) {
    mus.at(i) = stdVectorXd(dim, new_mu(i));
  }

  // set the value
  internal->setMu(mu_times, mus);

}

// set time/state varying mu
void TensorPhyloExternal::setMuTimeStateVarying(const VectorXd& new_mu_times, const MatrixXd& new_mu) {

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
  stdVectorXd mu_times = TensorPhyloUtils::EigenToStd(new_mu_times);

  // set the rate variable
  stdMatrixXd mus = TensorPhyloUtils::EigenToStd(new_mu);

  // set value
  internal->setMu(mu_times, mus);

}




/////////
// phi //
/////////

void TensorPhyloExternal::setPhiConstant(const double& new_phi) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
      stop("Error setting sampling rates. Rates must be strictly non-negative.");
    }

  }

  // set lambda times to empty
  stdVectorXd phi_times = stdVectorXd();

  // repeat the phis
  stdMatrixXd phis = stdMatrixXd(1, stdVectorXd(dim, new_phi));

  // set the value
  internal->setPhi(phi_times, phis);

}

// set state varying
void TensorPhyloExternal::setPhiStateVarying(const VectorXd& new_phi) {

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
  stdVectorXd phi_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_phi;

  stdMatrixXd phis = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setPhi(phi_times, phis);

}

// set time-varying phi
void TensorPhyloExternal::setPhiTimeVarying(const VectorXd& new_phi_times, const VectorXd& new_phi) {

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
  stdVectorXd phi_times = TensorPhyloUtils::EigenToStd(new_phi_times);

  // copy the time-varying rates per state
  stdMatrixXd phis(new_phi.size());
  for(size_t i = 0; i < new_phi.size(); ++i) {
    phis.at(i) = stdVectorXd(dim, new_phi(i));
  }

  // set the value
  internal->setPhi(phi_times, phis);

}

// set time/state varying phi
void TensorPhyloExternal::setPhiTimeStateVarying(const VectorXd& new_phi_times, const MatrixXd& new_phi) {

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
  stdVectorXd phi_times = TensorPhyloUtils::EigenToStd(new_phi_times);

  // set the rate variable
  stdMatrixXd phis = TensorPhyloUtils::EigenToStd(new_phi);

  // set value
  internal->setPhi(phi_times, phis);

}






///////////
// delta //
///////////

void TensorPhyloExternal::setDeltaConstant(const double& new_delta) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
      stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
    }

  }

  // set delta times to empty
  stdVectorXd delta_times = stdVectorXd();

  // repeat the deltas
  stdMatrixXd deltas = stdMatrixXd(1, stdVectorXd(dim, new_delta));

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set state varying
void TensorPhyloExternal::setDeltaStateVarying(const VectorXd& new_delta) {

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
  stdVectorXd delta_times = stdVectorXd();

  // create the rates
  MatrixXd ll(1, dim);
  ll.row(0) = new_delta;

  stdMatrixXd deltas = TensorPhyloUtils::EigenToStd(ll);

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set time-varying delta
void TensorPhyloExternal::setDeltaTimeVarying(const VectorXd& new_delta_times, const VectorXd& new_delta) {

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
  stdVectorXd delta_times = TensorPhyloUtils::EigenToStd(new_delta_times);

  // copy the time-varying rates per state
  stdMatrixXd deltas(new_delta.size());
  for(size_t i = 0; i < new_delta.size(); ++i) {
    deltas.at(i) = stdVectorXd(dim, new_delta(i));
  }

  // set the value
  internal->setDelta(delta_times, deltas);

}

// set time/state varying delta
void TensorPhyloExternal::setDeltaTimeStateVarying(const VectorXd& new_delta_times, const MatrixXd& new_delta) {

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
  stdVectorXd delta_times = TensorPhyloUtils::EigenToStd(new_delta_times);

  // set the rate variable
  stdMatrixXd deltas = TensorPhyloUtils::EigenToStd(new_delta);

  // set value
  internal->setDelta(delta_times, deltas);

}





/////////
// eta //
/////////

void TensorPhyloExternal::setEtaConstantEqual(const double& new_eta) {

  if ( safe ) {

    // check that it's valid rate
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta) == false ) {
      stop("Error setting transition rates. Rates must be strictly non-negative.");
    }

  }

  // set delta times to empty
  stdVectorXd eta_times = stdVectorXd();

  // create a matrix
  MatrixXd tmp       = MatrixXd::Constant(dim, dim, new_eta);
  tmp.diagonal()     = ((double)dim - 1) * VectorXd::Constant(dim, -new_eta);
  stdMatrixXd tmpStd = TensorPhyloUtils::EigenToStd(tmp);

  // create the vector of etas
  std::vector<stdMatrixXd> etas = std::vector<stdMatrixXd>(1, tmpStd);

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaConstantUnequal(const RateMatrix& new_eta) {

  if ( safe ) {

    // make sure this is a rate matrix
    if ( TensorPhyloUtils::isRateMatrix(new_eta) == false ) {
      stop("Error setting transition rates. Rate matrix must be square, and the rows must sum to zero.");
    }

    // check the dimensions
    if ( TensorPhyloUtils::hasDimensions(new_eta.getMatrix(), dim, dim) == false ) {
      stop("Error setting transition rates. Rate matrix must have one row and column per character state.");
    }

  }

  // set delta times to empty
  stdVectorXd eta_times = stdVectorXd();

  // create the matrix
  stdMatrixXd tmpStd = TensorPhyloUtils::EigenToStd( new_eta.getMatrix() );

  // create the vector of etas
  std::vector<stdMatrixXd> etas = std::vector<stdMatrixXd>(1, tmpStd);

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaTimeVaryingEqual(const VectorXd& new_eta_times, const VectorXd& new_eta) {

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
  stdVectorXd eta_times = TensorPhyloUtils::EigenToStd(new_eta_times);

  // make a rate matrix for each time
  std::vector<stdMatrixXd> etas(new_eta.size());
  for(size_t i = 0; i < new_eta.size(); ++i) {
    MatrixXd tmp   = MatrixXd::Constant(dim, dim, new_eta[i]);
    tmp.diagonal() = ((double)dim - 1) * VectorXd::Constant(dim, -new_eta[i]);
    etas.at(i) = TensorPhyloUtils::EigenToStd(tmp);
  }

  // set the value
  internal->setEta(eta_times, etas);

}

void TensorPhyloExternal::setEtaTimeVaryingUnequal(const VectorXd& new_eta_times, const RateMatrixList& new_eta) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta_times) == false ) {
      stop("Error setting transition rates. Times must be strictly non-negative.");
    }

    // loop over each matrix
    for(size_t i = 0; i < new_eta.size(); ++i) {

      // make sure this is a rate matrix
      if ( TensorPhyloUtils::isRateMatrix( new_eta.at(i) ) == false ) {
        stop("Error setting transition rates. Rate matrix must be square, and the rows must sum to zero.");
      }

      // check the dimensions
      if ( TensorPhyloUtils::hasDimensions(new_eta.at(i).getMatrix(), dim, dim) == false ) {
        stop("Error setting transition rates. Rate matrix must have one row and column per character state.");
      }

    }

  }

  // set the time variable
  stdVectorXd eta_times = TensorPhyloUtils::EigenToStd(new_eta_times);

  // create the "std" cube
  std::vector<stdMatrixXd> etas( new_eta.size() );
  for(size_t i = 0; i < new_eta.size(); ++i) {
    etas.at(i) = TensorPhyloUtils::EigenToStd(new_eta.at(i).getMatrix());
  }

  // set the value
  internal->setEta(eta_times, etas);

}



///////////
// omega //
///////////


void TensorPhyloExternal::setOmegaConstant(const CladoEvents& new_omega) {

  if ( safe ) {

    // check dimensions of matrix
    if ( TensorPhyloUtils::hasDimensions(new_omega, dim) == false ) {
      stop("Error setting cladogenetic events. Ancestral and daughter states must match the number of states for the data.");
    }

    // check that it's an event map
    if ( TensorPhyloUtils::isCladogeneticProbabilityMap(new_omega) == false ) {
      stop("Error setting cladogenetic events. Events for a given ancestral state must sum to 1.");
    }

  }

  // set omega times to empty
  stdVectorXd omega_times = stdVectorXd();

  // get the omegas
  eventMap_t omega( new_omega.getEvents() );
  std::vector<eventMap_t> omegas;
  omegas.push_back(omega);

  // set the value
  internal->setOmega(dim, omega_times, omegas);

}

void TensorPhyloExternal::setOmegaTimeVarying(const VectorXd& new_omega_times, const CladoEventsList& new_omegas) {

  if ( safe ) {

    // check that times are valid
    if ( TensorPhyloUtils::isStrictlyNonNegative(new_omega_times) == false ) {
      stop("Error setting cladogenetic events events. Times must be strictly non-negative.");
    }

    // loop over each element of list
    for(size_t i = 0; i < new_omegas.size(); ++i) {

      // check dimensions of matrix
      if ( TensorPhyloUtils::hasDimensions(new_omegas.at(i), dim) == false ) {
        stop("Error setting cladogenetic events. Ancestral and daughter states must match the number of states for the data.");
      }

      // check that it's an event map
      if ( TensorPhyloUtils::isCladogeneticProbabilityMap(new_omegas.at(i)) == false ) {
        stop("Error setting cladogenetic events. Events for a given ancestral state must sum to 1.");
      }

    }

  }

  // TODO: expose enums

  // set omega times to empty
  stdVectorXd omega_times = TensorPhyloUtils::EigenToStd(new_omega_times);

  // get the internal eventType_t from each omega
  // copy it, and emplace it in local omega
  std::vector<eventMap_t> omegas( new_omegas.size() );
  for(size_t i = 0; i < new_omegas.size(); ++i) {
    omegas.at(i) = eventMap_t( new_omegas.at(i).getEvents() );
  }

  // set the value
  internal->setOmega(dim, omega_times, omegas);

}

/////////////////////
// mass speciation //
/////////////////////

void TensorPhyloExternal::setUpsilonConstant(const VectorXd& new_upsilon_times, const VectorXd& new_upsilon) {

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
  stdVectorXd upsilon_times = TensorPhyloUtils::EigenToStd(new_upsilon_times);

  // copy the time-varying rates per state
  stdMatrixXd upsilons(new_upsilon.size());
  for(size_t i = 0; i < new_upsilon.size(); ++i) {
    upsilons.at(i) = stdVectorXd(dim, new_upsilon(i));
  }

  // set value
  internal->setMassSpeciationEvents(upsilon_times, upsilons);

}

void TensorPhyloExternal::setUpsilonStateVarying(const VectorXd& new_upsilon_times, const MatrixXd& new_upsilon) {

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
  stdVectorXd upsilon_times = TensorPhyloUtils::EigenToStd(new_upsilon_times);

  // set the rate variable
  stdMatrixXd upsilons = TensorPhyloUtils::EigenToStd(new_upsilon);

  // set value
  internal->setMassSpeciationEvents(upsilon_times, upsilons);

}

/////////////////////
// mass extinction //
/////////////////////

void TensorPhyloExternal::setGammaConstant(const VectorXd& new_gamma_times, const VectorXd& new_gamma) {

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
  stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

  // copy the time-varying rates per state
  stdMatrixXd gammas(new_gamma.size());
  for(size_t i = 0; i < new_gamma.size(); ++i) {
    gammas.at(i) = stdVectorXd(dim, new_gamma(i));
  }

  // set value
  internal->setMassSpeciationEvents(gamma_times, gammas);

}

void TensorPhyloExternal::setGammaStateVarying(const VectorXd& new_gamma_times, const MatrixXd& new_gamma) {

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
  stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

  // set the rate variable
  stdMatrixXd gammas = TensorPhyloUtils::EigenToStd(new_gamma);

  // set value
  internal->setMassSpeciationEvents(gamma_times, gammas);

}

// TODO: mass-extinction-induced state change



///////////////////
// mass sampling //
///////////////////

void TensorPhyloExternal::setRhoPresent(const double& new_rho) {

  if ( safe ) {

    // check that rho is a probability
    if ( TensorPhyloUtils::isProbability(new_rho) == false ) {
      stop("Error setting mass-sampling events. Magnitudes must be between 0 and 1 (inclusive).");
    }

  }

  // set time to the present
  stdVectorXd rho_times = stdVectorXd(1, 0);

  // create a matrix of sampling magnitudes
  MatrixXd tmp = MatrixXd::Constant(1, dim, new_rho);

  // set the matrix
  stdMatrixXd rhos = TensorPhyloUtils::EigenToStd(tmp);

  // set the value
  internal->setMassSamplingEvents(rho_times, rhos);

}

void TensorPhyloExternal::setRhoPresentStateVarying(const VectorXd& new_rho) {

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
  stdVectorXd rho_times = stdVectorXd(1, 0);

  // create the matrix
  stdVectorXd rho_vec = TensorPhyloUtils::EigenToStd(new_rho);
  stdMatrixXd rhos = stdMatrixXd(1, rho_vec);

  // set the value
  internal->setMassSamplingEvents(rho_times, rhos);

}

void TensorPhyloExternal::setRhoConstant(const VectorXd& new_rho_times, const VectorXd& new_rho) {

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
  stdVectorXd rho_times = TensorPhyloUtils::EigenToStd(new_rho_times);

  // copy the time-varying rates per state
  stdMatrixXd rhos(new_rho.size());
  for(size_t i = 0; i < new_rho.size(); ++i) {
    rhos.at(i) = stdVectorXd(dim, new_rho(i));
  }

  // set value
  internal->setMassSamplingEvents(rho_times, rhos);

}

void TensorPhyloExternal::setRhoStateVarying(const VectorXd& new_rho_times, const MatrixXd& new_rho) {

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
  stdVectorXd rho_times = TensorPhyloUtils::EigenToStd(new_rho_times);

  // set the rate variable
  stdMatrixXd rhos = TensorPhyloUtils::EigenToStd(new_rho);

  // set value
  internal->setMassSamplingEvents(rho_times, rhos);

}


///////////////////////////////
// mass destructive-sampling //
///////////////////////////////

void TensorPhyloExternal::setXiConstant(const VectorXd& new_xi_times, const VectorXd& new_xi) {

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
  stdVectorXd xi_times = TensorPhyloUtils::EigenToStd(new_xi_times);

  // copy the time-varying rates per state
  stdMatrixXd xis(new_xi.size());
  for(size_t i = 0; i < new_xi.size(); ++i) {
    xis.at(i) = stdVectorXd(dim, new_xi(i));
  }

  // set value
  internal->setMassDestrSamplingEvents(xi_times, xis);

}

void TensorPhyloExternal::setXiStateVarying(const VectorXd& new_xi_times, const MatrixXd& new_xi) {

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
  stdVectorXd xi_times = TensorPhyloUtils::EigenToStd(new_xi_times);

  // set the rate variable
  stdMatrixXd xis = TensorPhyloUtils::EigenToStd(new_xi);

  // set value
  internal->setMassDestrSamplingEvents(xi_times, xis);

}







// END: TensorPhyloExternal.cpp
