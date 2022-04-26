#ifndef TPEXT_H
#define TPEXT_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <boost/shared_ptr.hpp>
#include "Interface/DistributionHandlerImpl.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

RCPP_EXPOSED_CLASS(TensorPhyloExternal)
class TensorPhyloExternal {

  public:

    TensorPhyloExternal(size_t dim_);
    ~TensorPhyloExternal(){};

    // tree and data
    void setTree(const std::string &aNewickTree);
    void setData();

    // numerical settings
    void setSafeMode(bool safe_mode);
    void setNumberOfThreads(size_t nThreads);
    void setInitialDeltaT(double initDeltaT);
    void setLikelihoodApproximator(int approxVersion);
    void setSeed(size_t aSeed);

    // debugging
    void setDebugMode(int m);
    void setSyncMonitors(const std::vector< double > &synchMonitoring);
    void report() {
      Rcpp::Rcout << "A tensorphylo object at address " << internal << std::endl;
    }

    // model settings
    void setApplyTreeLikCorrection(bool doApply);
  	void setConditionalProbCompatibilityMode(bool setActive);
    void setConditionalProbabilityType(int condProb);

    // likelihood //
    double computeLogLikelihood();

    // root prior
    void     setRootPriorFlat();
    void     setRootPrior(VectorXd new_root_freq);

    // lambda
    void     setLambdaConstant(double new_lambda);
    void     setLambdaStateVarying(VectorXd new_lambda);
    void     setLambdaTimeVarying(VectorXd new_lambda_times, VectorXd new_lambda);
    void     setLambdaTimeStateVarying(VectorXd new_lambda_times, MatrixXd new_lambda);

    // mu
    void     setMuConstant(double new_mu);
    void     setMuStateVarying(VectorXd new_mu);
    void     setMuTimeVarying(VectorXd new_mu_times, VectorXd new_mu);
    void     setMuTimeStateVarying(VectorXd new_mu_times, MatrixXd new_mu);

    // phi
    void     setPhiConstant(double new_phi);
    void     setPhiStateVarying(VectorXd new_phi);
    void     setPhiTimeVarying(VectorXd new_phi_times, VectorXd new_phi);
    void     setPhiTimeStateVarying(VectorXd new_phi_times, MatrixXd new_phi);

    // delta
    void     setDeltaConstant(double new_delta);
    void     setDeltaStateVarying(VectorXd new_delta);
    void     setDeltaTimeVarying(VectorXd new_delta_times, VectorXd new_delta);
    void     setDeltaTimeStateVarying(VectorXd new_delta_times, MatrixXd new_delta);

    // eta
    void       setEtaConstantEqual(double new_eta);
    void       setEtaConstantUnequal(MatrixXd new_eta);
    void       setEtaTimeVaryingEqual(VectorXd new_eta_times, VectorXd new_eta);
    void       setEtaTimeVaryingUnequal(VectorXd new_eta_times, arma::cube new_eta);

    // // omega
    // // TODO: make omega somehow
    // void setOmega(size_t aNState, const NumericVector &times, const std::vector< eventMap_t > &omegas);

    // upsilon (mass speciation)
    void     setUpsilonConstant(VectorXd new_upsilon_times, VectorXd new_upsilon);
    void     setUpsilonStateVarying(VectorXd new_upsilon_times, MatrixXd new_upsilon);

    // gamma (mass-extinction events)
    void     setGammaConstant(VectorXd new_gamma_times, VectorXd new_gamma);
    void     setGammaStateVarying(VectorXd new_gamma_times, MatrixXd new_gamma);

    // zeta (mass-extinction-induced state change)

    // rho (mass-sampling events)
    void     setRhoPresent(double new_rho);
    void     setRhoPresentStateVarying(VectorXd new_rho);
    void     setRhoConstant(VectorXd new_rho_times, VectorXd new_rho);
    void     setRhoStateVarying(VectorXd new_rho_times, MatrixXd new_rho);

    // xi (mass destructive-sampling)
    void     setXiConstant(VectorXd new_xi_times, VectorXd new_xi);
    void     setXiStateVarying(VectorXd new_xi_times, MatrixXd new_xi);

  private:

    // pointer
    boost::shared_ptr<DistributionHandlerImpl> internal;

    // safe mode
    bool safe;

    // dimensions
    size_t dim;

};

#endif
