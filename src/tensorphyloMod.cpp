#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// using namespace Rcpp;
// using namespace Eigen;

///////////////////////
// tensorphylo class //
///////////////////////

#include "TensorPhyloExternal.h"

RCPP_MODULE(TensorPhyloMod) {

  // internal class
  class_<TensorPhyloExternal>("TensorPhylo")

    .constructor<int>()

    .method("setTree",      &TensorPhyloExternal::setTree)
    .method("setDebugMode", &TensorPhyloExternal::setDebugMode)
    .method("report",       &TensorPhyloExternal::report)

    // root prior rates
    .property("rootPrior",     &TensorPhyloExternal::getRootPrior, &TensorPhyloExternal::setRootPrior, "root prior")
    .method("updateRootPrior", &TensorPhyloExternal::updateRootPrior)

    // speciation rates
    .property("lambda",      &TensorPhyloExternal::getLambda,      &TensorPhyloExternal::setLambda,      "speciation rate")
    .property("lambdaTimes", &TensorPhyloExternal::getLambdaTimes, &TensorPhyloExternal::setLambdaTimes, "speciation rate change times")
    .method("updateLambdas", &TensorPhyloExternal::updateLambdas)

    // extinction rates
    .property("mu",      &TensorPhyloExternal::getMu,      &TensorPhyloExternal::setMu,      "extinction rate")
    .property("muTimes", &TensorPhyloExternal::getMuTimes, &TensorPhyloExternal::setMuTimes, "extinction rate change times")
    .method("updateMus", &TensorPhyloExternal::updateMus)

    // sampling rates
    .property("phi",      &TensorPhyloExternal::getPhi,      &TensorPhyloExternal::setPhi,      "sampling rate")
    .property("phiTimes", &TensorPhyloExternal::getPhiTimes, &TensorPhyloExternal::setPhiTimes, "sampling rate change times")
    .method("updatePhis", &TensorPhyloExternal::updatePhis)

    // destructive-sampling rates
    .property("delta",      &TensorPhyloExternal::getDelta,      &TensorPhyloExternal::setDelta,      "destructive-sampling rate")
    .property("deltaTimes", &TensorPhyloExternal::getDeltaTimes, &TensorPhyloExternal::setDeltaTimes, "destructive-sampling rate change times")
    .method("updateDeltas", &TensorPhyloExternal::updateDeltas)

    // anagenetic transition rates
    .property("eta",      &TensorPhyloExternal::getEta,      &TensorPhyloExternal::setEta,      "destructive-sampling rate")
    .property("etaTimes", &TensorPhyloExternal::getEtaTimes, &TensorPhyloExternal::setEtaTimes, "destructive-sampling rate change times")
    .method("updateEtas", &TensorPhyloExternal::updateEtas)

  ;

}










//
