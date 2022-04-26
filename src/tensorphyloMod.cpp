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

    // misc. (should be hidden)
    .method("setTree",      &TensorPhyloExternal::setTree)
    .method("setDebugMode", &TensorPhyloExternal::setDebugMode)
    .method("setSafeMode",  &TensorPhyloExternal::setSafeMode)
    .method("report",       &TensorPhyloExternal::report)

    // root frequency prior
    .method("getRootPrior",     &TensorPhyloExternal::getRootPrior)
    .method("setRootPriorFlat", &TensorPhyloExternal::setRootPriorFlat)
    .method("setRootPrior",     &TensorPhyloExternal::setRootPrior)

    // speciation rate
    .method("getLambda",                 &TensorPhyloExternal::getLambda)
    .method("getLambdaTimes",            &TensorPhyloExternal::getLambdaTimes)
    .method("setLambdaConstant",         &TensorPhyloExternal::setLambdaConstant)
    .method("setLambdaStateVarying",     &TensorPhyloExternal::setLambdaStateVarying)
    .method("setLambdaTimeVarying",      &TensorPhyloExternal::setLambdaTimeVarying)
    .method("setLambdaTimeStateVarying", &TensorPhyloExternal::setLambdaTimeStateVarying)

    // extinction rate
    .method("getMu",                 &TensorPhyloExternal::getMu)
    .method("getMuTimes",            &TensorPhyloExternal::getMuTimes)
    .method("setMuConstant",         &TensorPhyloExternal::setMuConstant)
    .method("setMuStateVarying",     &TensorPhyloExternal::setMuStateVarying)
    .method("setMuTimeVarying",      &TensorPhyloExternal::setMuTimeVarying)
    .method("setMuTimeStateVarying", &TensorPhyloExternal::setMuTimeStateVarying)

    // sampling rate
    .method("getPhi",                 &TensorPhyloExternal::getPhi)
    .method("getPhiTimes",            &TensorPhyloExternal::getPhiTimes)
    .method("setPhiConstant",         &TensorPhyloExternal::setPhiConstant)
    .method("setPhiStateVarying",     &TensorPhyloExternal::setPhiStateVarying)
    .method("setPhiTimeVarying",      &TensorPhyloExternal::setPhiTimeVarying)
    .method("setPhiTimeStateVarying", &TensorPhyloExternal::setPhiTimeStateVarying)

    // destructive-sampling rate
    .method("getDelta",                 &TensorPhyloExternal::getDelta)
    .method("getDeltaTimes",            &TensorPhyloExternal::getDeltaTimes)
    .method("setDeltaConstant",         &TensorPhyloExternal::setDeltaConstant)
    .method("setDeltaStateVarying",     &TensorPhyloExternal::setDeltaStateVarying)
    .method("setDeltaTimeVarying",      &TensorPhyloExternal::setDeltaTimeVarying)
    .method("setDeltaTimeStateVarying", &TensorPhyloExternal::setDeltaTimeStateVarying)

    // transition rates
    .method("getEta",                    &TensorPhyloExternal::getEta)
    .method("getEtaTimes",               &TensorPhyloExternal::getEtaTimes)
    .method("setEtaConstantEqual",       &TensorPhyloExternal::setEtaConstantEqual)
    .method("setEtaConstantUnequal",     &TensorPhyloExternal::setEtaConstantUnequal)
    .method("setEtaTimeVaryingEqual",    &TensorPhyloExternal::setEtaTimeVaryingEqual)
    .method("setEtaTimeVaryingUnequal",  &TensorPhyloExternal::setEtaTimeVaryingUnequal)

    // cladogenetic transition rates

    // mass-speciation events
    .method("getUpsilon",             &TensorPhyloExternal::getUpsilon)
    .method("getUpsilonTimes",        &TensorPhyloExternal::getUpsilonTimes)
    .method("setUpsilonConstant",     &TensorPhyloExternal::setUpsilonConstant)
    .method("setUpsilonStateVarying", &TensorPhyloExternal::setUpsilonStateVarying)

    // mass-extinction events
    .method("getGamma",             &TensorPhyloExternal::getGamma)
    .method("getGammaTimes",        &TensorPhyloExternal::getGammaTimes)
    .method("setGammaConstant",     &TensorPhyloExternal::setGammaConstant)
    .method("setGammaStateVarying", &TensorPhyloExternal::setGammaStateVarying)

    // mass-sampling events
    .method("getRho",                    &TensorPhyloExternal::getRho)
    .method("getRhoTimes",               &TensorPhyloExternal::getRhoTimes)
    .method("setRhoPresent",             &TensorPhyloExternal::setRhoPresent)
    .method("setRhoPresentStateVarying", &TensorPhyloExternal::setRhoPresentStateVarying)
    .method("setRhoConstant",            &TensorPhyloExternal::setRhoConstant)
    .method("setRhoStateVarying",        &TensorPhyloExternal::setRhoStateVarying)


    // mass-destructive-sampling events
    .method("getXi",             &TensorPhyloExternal::getXi)
    .method("getXiTimes",        &TensorPhyloExternal::getXiTimes)
    .method("setXiConstant",     &TensorPhyloExternal::setXiConstant)
    .method("setXiStateVarying", &TensorPhyloExternal::setXiStateVarying)


  ;

}










//
