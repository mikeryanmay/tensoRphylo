#ifndef TPEXT_H
#define TPEXT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/shared_ptr.hpp>
#include "Interface/DistributionHandlerImpl.h"
#include "CladoEvents.h"
#include "RateMatrix.h"
#include "ProbabilityMatrix.h"

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
    void setRootPriorFlat();
    void setRootPrior(const VectorXd& new_root_freq);

    // lambda
    void setLambdaConstant(const double& new_lambda);
    void setLambdaStateVarying(const VectorXd& new_lambda);
    void setLambdaTimeVarying(const VectorXd& new_lambda_times, const VectorXd& new_lambda);
    void setLambdaTimeStateVarying(const VectorXd& new_lambda_times, const MatrixXd& new_lambda);

    // mu
    void setMuConstant(const double& new_mu);
    void setMuStateVarying(const VectorXd& new_mu);
    void setMuTimeVarying(const VectorXd& new_mu_times, const VectorXd& new_mu);
    void setMuTimeStateVarying(const VectorXd& new_mu_times, const MatrixXd& new_mu);

    // phi
    void setPhiConstant(const double& new_phi);
    void setPhiStateVarying(const VectorXd& new_phi);
    void setPhiTimeVarying(const VectorXd& new_phi_times, const VectorXd& new_phi);
    void setPhiTimeStateVarying(const VectorXd& new_phi_times, const MatrixXd& new_phi);

    // delta
    void setDeltaConstant(const double& new_delta);
    void setDeltaStateVarying(const VectorXd& new_delta);
    void setDeltaTimeVarying(const VectorXd& new_delta_times, const VectorXd& new_delta);
    void setDeltaTimeStateVarying(const VectorXd& new_delta_times, const MatrixXd& new_delta);

    // eta
    void setEtaConstantEqual(const double& new_eta);
    void setEtaConstantUnequal(const RateMatrix& new_eta);
    void setEtaTimeVaryingEqual(const VectorXd& new_eta_times, const VectorXd& new_eta);
    void setEtaTimeVaryingUnequal(const VectorXd& new_eta_times, const RateMatrixList& new_eta);

    // omega
    void setOmegaConstant(const CladoEvents& new_omega);
    void setOmegaTimeVarying(const VectorXd& new_omega_times, const CladoEventsList& new_omegas);

    // upsilon (mass speciation)
    void setUpsilonConstant(const VectorXd& new_upsilon_times, const VectorXd& new_upsilon);
    void setUpsilonStateVarying(const VectorXd& new_upsilon_times, const MatrixXd& new_upsilon);

    // gamma (mass-extinction events)
    void setGammaConstant(const VectorXd& new_gamma_times, const VectorXd& new_gamma);
    void setGammaStateVarying(const VectorXd& new_gamma_times, const MatrixXd& new_gamma);

    // zeta (mass-extinction-induced state change)
    void setGammaAndZetaConstant(const VectorXd& new_gamma_times, const VectorXd& new_gamma, const ProbabilityMatrixList& new_zeta);
    void setGammaAndZetaStateVarying(const VectorXd& new_gamma_times, const MatrixXd& new_gamma, const ProbabilityMatrixList& new_zeta);

    // zeta (without mass-extinction events)
    void setZeta(const VectorXd& new_gamma_times, const ProbabilityMatrixList& new_zeta);

    // rho (mass-sampling events)
    void setRhoPresent(const double& new_rho);
    void setRhoPresentStateVarying(const VectorXd& new_rho);
    void setRhoConstant(const VectorXd& new_rho_times, const VectorXd& new_rho);
    void setRhoStateVarying(const VectorXd& new_rho_times, const MatrixXd& new_rho);

    // xi (mass destructive-sampling)
    void setXiConstant(const VectorXd& new_xi_times, const VectorXd& new_xi);
    void setXiStateVarying(const VectorXd& new_xi_times, const MatrixXd& new_xi);

  private:

    // pointer
    boost::shared_ptr<DistributionHandlerImpl> internal;

    // safe mode
    bool safe;

    // dimensions
    size_t dim;

};

#endif
