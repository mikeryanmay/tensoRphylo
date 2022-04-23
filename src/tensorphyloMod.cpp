#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/shared_ptr.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif
using namespace Rcpp;

///////////////////////
// tensorphylo class //
///////////////////////

#include "Interface/DistributionHandler.h"
#include "Interface/DistributionHandlerImpl.h"

// [[Rcpp::export]]
void testTPDistHandler() {

  using namespace TensorPhylo::Interface;
  boost::shared_ptr<DistributionHandlerImpl> tmp = DistributionHandlerImpl::create();

  return;

}


// RCPP_MODULE(tensorphyloMod) {
//
//   // the external object
//   class_<External>("External")
//     .constructor<int>()
//     .method("observe",          &External::observe)
//     .method("updateDependents", &External::updateDependents)
//     .method("report",           &External::report)
//   ;
//
//   // generic functions
//   function("printInternal", &printInternal, "");
//
//
// }
