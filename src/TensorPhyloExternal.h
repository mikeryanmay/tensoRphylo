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

    TensorPhyloExternal();
    TensorPhyloExternal(size_t dim_);
    ~TensorPhyloExternal(){};

    // tree and data
    void setTree(const std::string &aNewickTree);
    void setData();

    // numerical settings
    void setNumberOfThreads(size_t nThreads);
    void setInitialDeltaT(double initDeltaT);
    void setLikelihoodApproximator(int approxVersion);

    // debugging
    void setDebugMode(int m);
    void setSyncMonitors(const std::vector< double > &synchMonitoring);
    void setSeed(size_t aSeed);
    void report() {
      Rcpp::Rcout << "A tensorphylo object at address " << internal << std::endl;
    }

    // model settings
    void setApplyTreeLikCorrection(bool doApply);
  	void setConditionalProbCompatibilityMode(bool setActive);
    void setConditionalProbabilityType(int condProb);

    // root prior
    VectorXd getRootPrior();
    void     setRootPrior(VectorXd new_root_freq);
    void     updateRootPrior();

    // lambda
    MatrixXd getLambda();
    void     setLambda(MatrixXd new_lambda);
    VectorXd getLambdaTimes();
    void     setLambdaTimes(VectorXd new_lambda);
    void     updateLambdas();

    // mu
    MatrixXd getMu();
    void     setMu(MatrixXd new_mu);
    VectorXd getMuTimes();
    void     setMuTimes(VectorXd new_mu);
    void     updateMus();

    // phi
    MatrixXd getPhi();
    void     setPhi(MatrixXd new_phi);
    VectorXd getPhiTimes();
    void     setPhiTimes(VectorXd new_phi);
    void     updatePhis();

    // delta
    MatrixXd getDelta();
    void     setDelta(MatrixXd new_delta);
    VectorXd getDeltaTimes();
    void     setDeltaTimes(VectorXd new_delta);
    void     updateDeltas();

    // eta
    arma::cube getEta();
    void       setEta(arma::cube new_eta);
    VectorXd   getEtaTimes();
    void       setEtaTimes(VectorXd new_eta);
    void       updateEtas();

    // omega
    // TODO: make omega somehow
    void setOmega(size_t aNState, const NumericVector &times, const std::vector< eventMap_t > &omegas);

    // mass speciation
    void setMassSpeciationEvents(const NumericVector &massSpeciationTimes, const NumericMatrix &massSpeciationProb);

    // mass extinction
    void setMassExtinctionEvents(const NumericVector &massExtinctionTimes, const NumericMatrix &massExtinctionProb);

    // mass-extinction-induced state change
    void setMassExtinctionStateChangeProb(const std::vector< NumericMatrix> &massExtinctionStateChangeProb);

    // mass-sampling events
    void setMassSamplingEvents(const NumericVector &massSamplingTimes, const NumericMatrix &massSamplingProb);

    // mass-destructive-sampling events
    void setMassDestrSamplingEvents(const NumericVector &massDestrSamplingTimes, const NumericMatrix &massDestrSamplingProb);

    // likelihoods
    double computeLogLikelihood();

    // TODO: history stuff

  private:

    // dimensions
    size_t dim;

    // converters
    stdVectorXd              EigenToStd(VectorXd eig_vec);
    stdMatrixXd              EigenToStd(MatrixXd eig_mat);
    std::vector<stdMatrixXd> ArmaToStd(arma::cube arma_cube);

    // parameters
    bool     rf_dirty;
    VectorXd rootFrequency;

    bool     lambda_dirty;
    MatrixXd lambdas;
    VectorXd lambda_times;

    bool     mu_dirty;
    MatrixXd mus;
    VectorXd mu_times;

    bool     phi_dirty;
    MatrixXd phis;
    VectorXd phi_times;

    bool     delta_dirty;
    MatrixXd deltas;
    VectorXd delta_times;

    bool       eta_dirty;
    arma::cube etas;
    VectorXd   eta_times;

    // pointer
    boost::shared_ptr<DistributionHandlerImpl> internal;



};

#endif
