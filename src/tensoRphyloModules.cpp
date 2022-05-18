#include <Rcpp.h>
#include <RcppEigen.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

#include "TensorPhyloExternal.h"
#include "DataReader.h"

////////////////////
// expose classes //
////////////////////

RCPP_EXPOSED_CLASS(CladoEvents)
RCPP_EXPOSED_CLASS(CladoEventsList)
RCPP_EXPOSED_CLASS(ProbabilityMatrix)
RCPP_EXPOSED_CLASS(ProbabilityMatrixList)
RCPP_EXPOSED_CLASS(RateMatrix)
RCPP_EXPOSED_CLASS(RateMatrixList)
RCPP_EXPOSED_CLASS(TensorPhyloExternal)

///////////////////////
// data reader class //
///////////////////////

using namespace Rcpp;

RCPP_MODULE(DataReaderMod) {

  Rcpp::function("readDelimitedData", &DataReader::readDelimitedData);
  Rcpp::function("readNexusData",     &DataReader::readNexusData);

}

///////////////////////////////////
// create the R interface module //
///////////////////////////////////

RCPP_MODULE(TensorPhyloMod) {

  // RateMatrix class
  class_<RateMatrix>("RateMatrix")
    .constructor<size_t, double>()
    .constructor<MatrixXd>()
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
    .constructor<size_t, double>()
    .constructor<MatrixXd>()
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
  class_<TensorPhyloExternal>("TensorPhyloInstance")

    .constructor<size_t>()
    .constructor<List, std::string, NumericMatrix>()

    // data and tree
    .method("setTree", &TensorPhyloExternal::setTree)
    .method("setData", &TensorPhyloExternal::setData)

    // the main event
    .method("computeLogLikelihood",          &TensorPhyloExternal::computeLogLikelihood)
    .method("drawStochasticMaps",            &TensorPhyloExternal::drawStochasticMaps)
    .method("drawBranchRates",               &TensorPhyloExternal::drawBranchRates)

    // misc.
    .method("setNumberOfThreads",              &TensorPhyloExternal::setNumberOfThreads)
    .method("setApplyTreeLikCorrection",       &TensorPhyloExternal::setApplyTreeLikCorrection)
    .method("setQuasistationaryFrequencyMode", &TensorPhyloExternal::setQuasistationaryFrequencyMode)
    .method("setConditionalProbabilityType",   &TensorPhyloExternal::setConditionalProbabilityType)
    .method("setLikelihoodApproximator",       &TensorPhyloExternal::setLikelihoodApproximator)
    .method("setIntegrationScheme",            &TensorPhyloExternal::setIntegrationScheme)
    .method("setDebugMode",                    &TensorPhyloExternal::setDebugMode)
    .method("setSafeMode",                     &TensorPhyloExternal::setSafeMode)
    .method("show",                            &TensorPhyloExternal::report)

    // root frequency prior
    .method("setRootPriorFlat",            &TensorPhyloExternal::setRootPriorFlat)
    .method("setRootPrior",                &TensorPhyloExternal::setRootPrior)
    .method("getQuasiStationaryFrequency", &TensorPhyloExternal::getQuasiStationaryFrequency)

    // speciation rate
    .method("setLambdaConstant",           &TensorPhyloExternal::setLambdaConstant)
    .method("setLambdaStateDependent",     &TensorPhyloExternal::setLambdaStateDependent)
    .method("setLambdaTimeDependent",      &TensorPhyloExternal::setLambdaTimeDependent)
    .method("setLambdaTimeStateDependent", &TensorPhyloExternal::setLambdaTimeStateDependent)

    // extinction rate
    .method("setMuConstant",           &TensorPhyloExternal::setMuConstant)
    .method("setMuStateDependent",     &TensorPhyloExternal::setMuStateDependent)
    .method("setMuTimeDependent",      &TensorPhyloExternal::setMuTimeDependent)
    .method("setMuTimeStateDependent", &TensorPhyloExternal::setMuTimeStateDependent)

    // sampling rate
    .method("setPhiConstant",           &TensorPhyloExternal::setPhiConstant)
    .method("setPhiStateDependent",     &TensorPhyloExternal::setPhiStateDependent)
    .method("setPhiTimeDependent",      &TensorPhyloExternal::setPhiTimeDependent)
    .method("setPhiTimeStateDependent", &TensorPhyloExternal::setPhiTimeStateDependent)

    // destructive-sampling rate
    .method("setDeltaConstant",           &TensorPhyloExternal::setDeltaConstant)
    .method("setDeltaStateDependent",     &TensorPhyloExternal::setDeltaStateDependent)
    .method("setDeltaTimeDependent",      &TensorPhyloExternal::setDeltaTimeDependent)
    .method("setDeltaTimeStateDependent", &TensorPhyloExternal::setDeltaTimeStateDependent)

    // transition rates
    .method("setEtaConstantEqual",         &TensorPhyloExternal::setEtaConstantEqual)
    .method("setEtaConstantUnequal",       &TensorPhyloExternal::setEtaConstantUnequal)
    .method("setEtaTimeDependentEqual",    &TensorPhyloExternal::setEtaTimeDependentEqual)
    .method("setEtaTimeDependentUnequal",  &TensorPhyloExternal::setEtaTimeDependentUnequal)

    // cladogenetic transition rates
    .method("setOmegaConstant",      &TensorPhyloExternal::setOmegaConstant)
    .method("setOmegaTimeDependent", &TensorPhyloExternal::setOmegaTimeDependent)

    // mass-speciation events
    .method("setUpsilonConstant",       &TensorPhyloExternal::setUpsilonConstant)
    .method("setUpsilonStateDependent", &TensorPhyloExternal::setUpsilonStateDependent)

    // mass-extinction events
    .method("setGammaConstant",            &TensorPhyloExternal::setGammaConstant)
    .method("setGammaStateDependent",        &TensorPhyloExternal::setGammaStateDependent)
    .method("setGammaAndZetaConstant",     &TensorPhyloExternal::setGammaAndZetaConstant)
    .method("setGammaAndZetaStateDependent", &TensorPhyloExternal::setGammaAndZetaStateDependent)
    .method("setZeta",                     &TensorPhyloExternal::setZeta)

    // mass-sampling events
    .method("setRhoPresent",               &TensorPhyloExternal::setRhoPresent)
    .method("setRhoPresentStateDependent", &TensorPhyloExternal::setRhoPresentStateDependent)
    .method("setRhoConstant",              &TensorPhyloExternal::setRhoConstant)
    .method("setRhoStateDependent",        &TensorPhyloExternal::setRhoStateDependent)

    // mass-destructive-sampling events
    .method("setXiConstant",       &TensorPhyloExternal::setXiConstant)
    .method("setXiStateDependent", &TensorPhyloExternal::setXiStateDependent)

  ;

}










//
