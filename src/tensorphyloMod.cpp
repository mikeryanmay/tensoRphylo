#include <Rcpp.h>
#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

///////////////////////
// tensorphylo class //
///////////////////////

#include "CladoEvents.h"
#include "RateMatrix.h"
#include "ProbabilityMatrix.h"
#include "TensorPhyloExternal.h"

///////////////////////////////////
// create the R interface module //
///////////////////////////////////

RCPP_MODULE(TensorPhyloMod) {

  // RateMatrix class
  class_<RateMatrix>("RateMatrix")
    .constructor<size_t>()
    .constructor<size_t, double>()
    .method("getMatrix",             &RateMatrix::getMatrix)
    .method("getRateOneIndexed",     &RateMatrix::getRate)
    .method("setRateOneIndexed",     &RateMatrix::setRate)
    .method("show",                  &RateMatrix::show)
  ;

  // RateMatrix list class
  class_<RateMatrixList>("RateMatrixList")
    .constructor<size_t>()
    .method("addRateMatrix", &RateMatrixList::addRateMatrix)
    .method("show",          &RateMatrixList::show)
  ;

  // ProbabilityMatrix class
  class_<ProbabilityMatrix>("ProbabilityMatrix")
    .constructor<size_t>()
    .method("getMatrix",                 &ProbabilityMatrix::getMatrix)
    .method("getProbabilityOneIndexed",  &ProbabilityMatrix::getProbability)
    .method("setProbabilityOneIndexed",  &ProbabilityMatrix::setProbability)
    .method("show",                      &ProbabilityMatrix::show)
  ;

  // ProbabilityMatrix list class
  class_<ProbabilityMatrixList>("ProbabilityMatrixList")
    .constructor<size_t>()
    .method("addProbabilityMatrix", &ProbabilityMatrixList::addProbabilityMatrix)
    .method("show",                 &ProbabilityMatrixList::show)
  ;

  // CladoEvents class
  class_<CladoEvents>("CladoEvents")
    .constructor<size_t>()
    .method("getEventOneIndexed", &CladoEvents::getEvent)
    .method("setEventOneIndexed", &CladoEvents::setEvent)
    .method("show",               &CladoEvents::show)
  ;

  // CladoEvents list class
  class_<CladoEventsList>("CladoEventsList")
    .constructor<size_t>()
    .method("addCladoEvents", &CladoEventsList::addCladoEvents)
    .method("show",           &CladoEventsList::show)
  ;

  // internal class
  class_<TensorPhyloExternal>("TensorPhylo")

    .constructor<size_t>()
    .constructor<List, std::string, NumericMatrix>()

    // the main event
    .method("computeLogLikelihood", &TensorPhyloExternal::computeLogLikelihood)
    .method("drawStochasticMap",    &TensorPhyloExternal::drawStochasticMap)

    // misc.
    .method("setDebugMode", &TensorPhyloExternal::setDebugMode)
    .method("setSafeMode",  &TensorPhyloExternal::setSafeMode)
    .method("show",         &TensorPhyloExternal::report)

    // root frequency prior
    .method("setRootPriorFlat", &TensorPhyloExternal::setRootPriorFlat)
    .method("setRootPrior",     &TensorPhyloExternal::setRootPrior)

    // speciation rate
    .method("setLambdaConstant",         &TensorPhyloExternal::setLambdaConstant)
    .method("setLambdaStateVarying",     &TensorPhyloExternal::setLambdaStateVarying)
    .method("setLambdaTimeVarying",      &TensorPhyloExternal::setLambdaTimeVarying)
    .method("setLambdaTimeStateVarying", &TensorPhyloExternal::setLambdaTimeStateVarying)

    // extinction rate
    .method("setMuConstant",         &TensorPhyloExternal::setMuConstant)
    .method("setMuStateVarying",     &TensorPhyloExternal::setMuStateVarying)
    .method("setMuTimeVarying",      &TensorPhyloExternal::setMuTimeVarying)
    .method("setMuTimeStateVarying", &TensorPhyloExternal::setMuTimeStateVarying)

    // sampling rate
    .method("setPhiConstant",         &TensorPhyloExternal::setPhiConstant)
    .method("setPhiStateVarying",     &TensorPhyloExternal::setPhiStateVarying)
    .method("setPhiTimeVarying",      &TensorPhyloExternal::setPhiTimeVarying)
    .method("setPhiTimeStateVarying", &TensorPhyloExternal::setPhiTimeStateVarying)

    // destructive-sampling rate
    .method("setDeltaConstant",         &TensorPhyloExternal::setDeltaConstant)
    .method("setDeltaStateVarying",     &TensorPhyloExternal::setDeltaStateVarying)
    .method("setDeltaTimeVarying",      &TensorPhyloExternal::setDeltaTimeVarying)
    .method("setDeltaTimeStateVarying", &TensorPhyloExternal::setDeltaTimeStateVarying)

    // transition rates
    .method("setEtaConstantEqual",       &TensorPhyloExternal::setEtaConstantEqual)
    .method("setEtaConstantUnequal",     &TensorPhyloExternal::setEtaConstantUnequal)
    .method("setEtaTimeVaryingEqual",    &TensorPhyloExternal::setEtaTimeVaryingEqual)
    .method("setEtaTimeVaryingUnequal",  &TensorPhyloExternal::setEtaTimeVaryingUnequal)

    // cladogenetic transition rates
    .method("setOmegaConstant",    &TensorPhyloExternal::setOmegaConstant)
    .method("setOmegaTimeVarying", &TensorPhyloExternal::setOmegaTimeVarying)

    // mass-speciation events
    .method("setUpsilonConstant",     &TensorPhyloExternal::setUpsilonConstant)
    .method("setUpsilonStateVarying", &TensorPhyloExternal::setUpsilonStateVarying)

    // mass-extinction events
    .method("setGammaConstant",            &TensorPhyloExternal::setGammaConstant)
    .method("setGammaStateVarying",        &TensorPhyloExternal::setGammaStateVarying)
    .method("setGammaAndZetaConstant",     &TensorPhyloExternal::setGammaAndZetaConstant)
    .method("setGammaAndZetaStateVarying", &TensorPhyloExternal::setGammaAndZetaStateVarying)
    .method("setZeta",                     &TensorPhyloExternal::setZeta)

    // mass-sampling events
    .method("setRhoPresent",             &TensorPhyloExternal::setRhoPresent)
    .method("setRhoPresentStateVarying", &TensorPhyloExternal::setRhoPresentStateVarying)
    .method("setRhoConstant",            &TensorPhyloExternal::setRhoConstant)
    .method("setRhoStateVarying",        &TensorPhyloExternal::setRhoStateVarying)

    // mass-destructive-sampling events
    .method("setXiConstant",     &TensorPhyloExternal::setXiConstant)
    .method("setXiStateVarying", &TensorPhyloExternal::setXiStateVarying)


  ;

}










//
