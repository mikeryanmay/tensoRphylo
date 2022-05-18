#ifndef TPEXT_H
#define TPEXT_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/shared_ptr.hpp>

// include openMP, if available
#ifdef _OPENMP
  #include <omp.h>
#endif

// now for tensorphylo stuff
#include "Interface/DistributionHandlerImpl.h"
#include "TensorPhyloUtils.h"
#include "CladoEvents.h"
#include "RateMatrix.h"
#include "ProbabilityMatrix.h"

using namespace Rcpp;
using namespace Eigen;
using namespace TensorPhylo::Interface;

class TensorPhyloExternal {

  public:

    TensorPhyloExternal(size_t dim_) :
      ready(false),
      safe(true),
      dim(dim_),
      has_omega(false) {

      // make some state labels
      state_labels.resize(dim);
      for(size_t i = 0; i < dim; ++i) {
        state_labels.at(i) = std::to_string(i);
      }

      // initialize
      init();

    };

    TensorPhyloExternal(List phylo, std::string newick, NumericMatrix data) :
      ready(false),
      safe(true),
      has_omega(false) {

      // get the dimensions
      dim = data.ncol();

      // make some state labels
      state_labels.resize(dim);
      for(size_t i = 0; i < dim; ++i) {
        state_labels.at(i) = std::to_string(i);
      }

      // initialize
      init();

      // store some stuff
      phy = phylo;

      // initialize the tree
      setTree(newick);

      // initialize data
      setData(data);

      // the distribution is ready
      ready = true;

    }

    ///////////////////
    // tree and data //
    ///////////////////

    void setTree(const std::string &aNewickTree) {
      internal->setTree(aNewickTree);
    }

    void setData(const NumericMatrix& aProbMatrix) {

      // get the state labels
      state_labels = as<std::vector<std::string>>(colnames(aProbMatrix));

      // get the taxa
      CharacterVector taxa = rownames(aProbMatrix);

      // create the map
      std::vector<std::string> labels;
      std::map<std::string, std::vector<double>> probabilityMap;

      // loop over each taxon
      for(size_t i = 0; i < aProbMatrix.nrow(); ++i) {

        // the first element is the taxon name
        std::string taxon = as<std::string>(taxa[i]);

        // the second element is the vector of probabilities
        NumericVector nv = aProbMatrix.row(i);
        std::vector<double> prob = as<std::vector<double> >( nv );

        // insert the data
        labels.push_back(taxon);
        probabilityMap.emplace(taxon, prob);

      }

      // send the data to the internal object
      internal->setData(labels, probabilityMap);

    }

    ////////////////////////
    // numerical settings //
    ////////////////////////

    void setSafeMode(bool safe_mode) {
      safe = safe_mode;
    }

    void setNumberOfThreads(size_t nThreads) {
      #ifdef _OPENMP
        internal->setNumberOfThreads(nThreads);
      #else
        stop("OpenMP not found.");
      #endif
    }

    void setInitialDeltaT(double initDeltaT) {
      internal->setInitialDeltaT(initDeltaT);
    }

    void setLikelihoodApproximator(int approxVersion) {
      internal->setLikelihoodApproximator( (approximatorVersion_t)approxVersion );
    }

    void setIntegrationScheme(int aIntScheme) {
      internal->setIntegrationScheme( (integrationScheme_t)aIntScheme );
    }

    void setSeed(size_t aSeed) {
      internal->setSeed(aSeed);
    }

    ///////////////
    // debugging //
    ///////////////

    void setDebugMode(int m) {
      internal->setDebugMode( (debugMode_t)m );
    }

    void setSyncMonitors(const std::vector< double > &synchMonitoring) {
      internal->setSyncMonitors(synchMonitoring);
    }

    void report() {
      Rcout << "A TensorPhyloInstance object at address <" << internal << ">" << std::endl;
    }

    ////////////////////
    // model settings //
    ////////////////////

    void setApplyTreeLikCorrection(bool doApply) {
      internal->setApplyTreeLikCorrection(doApply);
    }

    void setQuasistationaryFrequencyMode(bool setActive) {
      internal->setQuasistationaryFrequencyMode(setActive);
    }

    void setConditionalProbCompatibilityMode(bool setActive) {
      internal->setConditionalProbCompatibilityMode(setActive);
    }

    void setConditionalProbabilityType(int condProb) {
      internal->setConditionalProbabilityType( (conditionalProbability_t)condProb );
    }

    ////////////////
    // likelihood //
    ////////////////

    double computeLogLikelihood() {
      if ( ready == false ) {
        stop("Distribution not ready. Needs data and parameter.");
      }
      return internal->computeLogLikelihood();
    }

    ////////////////////////
    // stochastic mapping //
    ////////////////////////

    List drawStochasticMaps(size_t reps) {

      if ( ready == false ) {
        stop("Distribution not ready. Needs data and parameter.");
      }

      // create a list for return
      List res(reps);

      // loop over reps
      Rcout << "Simulating stochastic maps. Please be patient." << std::endl;
      size_t i = 0; // replicate incrementer
      try {

        for(; i < reps; ++i) {

          // initialize the phylo
          List this_phy = phy;

          // get a stochastic map
          // mapHistories_t map = maps.at(i);
          mapHistories_t map = internal->drawHistory();

          // create the iterator for branch histories
          mapHistories_t::iterator it = map.begin();

          // manually increment the iterator to ignore the stem
          it++;

          // make containers for the branch histories and the mapped.edge
          List histories(map.size() - 1);
          NumericMatrix dwell_time(map.size() - 1, dim);
          colnames(dwell_time) = wrap(state_labels);

          // loop over the remaining branch histories
          size_t history_index = 0;
          for(; it != map.end(); ++it) {

            // get the history for the branch
            mapHistoriesVal_t this_history = it->second;

            // this is a vector of pairs <time, state>
            std::vector<double>      branch_times;
            std::vector<std::string> branch_states;
            for(mapHistoriesVal_t::iterator jt = this_history.begin(); jt != this_history.end(); ++jt) {

              // add the branch duration with a name
              branch_times.push_back( jt->first );
              branch_states.push_back( state_labels.at(jt->second) );

              // increment the dwell time
              dwell_time(history_index, jt->second) += jt->first;

            } // end loop over branch events

            // create a named vector
            NumericVector char_history = wrap(branch_times);
            char_history.names() = wrap(branch_states);

            // insert the history
            histories[history_index++] = char_history;

          } // end loop over branch histories

          // attach the maps to the tree
          this_phy["maps"] = histories;
          this_phy["mapped.edge"] = dwell_time;
          this_phy.attr("class") = CharacterVector::create("simmap","phylo");

          // store the tree
          res.at(i) = this_phy;

          // check for user interrupt
          checkUserInterrupt();

          // increment the progress bar
          // double percent = 100.0 * (double)i / (double)reps;
          double percent = 100.0 * (double)(i + 1) / (double)reps;
          int perc = floor(percent);
          int frac = floor(0.4 * percent);
          Rprintf("\r");
          Rcout << "[";
          for(size_t j = 0; j < 40; ++j) {
            if ( j <= frac ) {
              Rcout << "=";
            } else {
              Rcout << " ";
            }
          }
          Rcout << "] ";
          Rprintf("%d%%", perc);
          Rprintf("\r");

          // only flush the console if we're not on windows
          #if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
          R_FlushConsole();
          #endif

        } // end loop over replicates

      } catch ( Rcpp::internal::InterruptedException& e ) {

        // terminate the progress bar
        Rcout << std::endl;

        // subset the results so far
        res = res[Range(0, i)];

        // make sure the list is a multisimmap
        res.attr("class") = CharacterVector::create("multiSimmap","multiPhylo");

        // return the simulations
        return res;

      } // catch user interruption

      // terminate the progress bar
      Rcout << std::endl;

      // make sure the list is a multisimmap
      res.attr("class") = CharacterVector::create("multiSimmap","multiPhylo");

      // return the simulations
      return res;

    }

    List drawBranchRates(size_t reps) {

      if ( ready == false ) {
        stop("Distribution not ready. Needs data and parameter.");
      }

      // initialize the return list
      List res(reps);

      NumericVector edges = phy["edge.length"];
      size_t nrow = edges.size();

      Rcout << "Simulating branch rates. Please be patient." << std::endl;
      size_t i = 0; // replicate incrementer
      try {

        // loop over simulation
        for(; i < reps; ++i) {

          // create the containers
          std::vector<double> lambda;
          std::vector<double> mu;
          std::vector<double> phi;
          std::vector<double> delta;
          std::vector<long>   num;

          // do a simulation
          internal->drawHistoryAndComputeRates(lambda, mu, phi, delta, num);

          // translate to Rcpp objects
          NumericVector l = wrap(lambda);
          NumericVector m = wrap(mu);
          NumericVector p = wrap(phi);
          NumericVector d = wrap(delta);
          NumericVector n = wrap(num);

          // create the container
          NumericMatrix rates(l.size(), 5);
          rates.column(0) = l;
          rates.column(1) = m;
          rates.column(2) = p;
          rates.column(3) = d;
          rates.column(4) = n;

          // store the result
          res[i] = rates;

          // check for user interrupt
          checkUserInterrupt();

          // increment the progress bar
          // double percent;
          // if ( i == reps ) {
          //   percent = 100.0;
          // } else {
          //   percent = 100.0 * (double)i / (double)reps;
          // }
          double percent = 100.0 * (double)(i + 1) / (double)reps;
          int perc = floor(percent);
          int frac = floor(0.4 * percent);
          Rprintf("\r");
          Rcout << "[";
          for(size_t j = 0; j <= 40; ++j) {
            if ( j <= frac ) {
              Rcout << "=";
            } else {
              Rcout << " ";
            }
          }
          Rcout << "] ";
          Rprintf("%d%%", perc);
          Rprintf("\r");

          // only flush the console if we're not on windows
          #if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
          R_FlushConsole();
          #endif

        } // end loop over replicates

        // terminate the progress bar
        Rcout << std::endl;

      } catch ( Rcpp::internal::InterruptedException& e ) {

        // terminate early

        // terminate the progress bar
        Rcout << std::endl;

        // subset the results so far
        res = res[Range(0, i)];

        // return the simulations
        return res;

      }

      // return the matrices
      return res;

    }

    ////////////////
    // root prior //
    ////////////////

    void setRootPriorFlat() {

      // make a flat vector
      stdVectorXd root_frequency = stdVectorXd(dim, 1.0 / (double)dim);

      // update
      internal->setRootPrior(root_frequency);

    }

    void setRootPrior(const VectorXd& new_root_freq) {

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

    void setLambdaConstant(const double& new_lambda) {

      if ( safe ) {

        // check that time is valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_lambda) == false ) {
          stop("Error setting speciation rates. Rates must be strictly non-negative.");
        }

      }

      // set lambda times to empty
      stdVectorXd lambda_times = stdVectorXd();

      // repeat the lambdas
      stdMatrixXd lambdas = TensorPhyloUtils::ScalarToStdMat(new_lambda, dim);

      // set the value
      internal->setLambda(lambda_times, lambdas);

    }

    // set state dependent
    void setLambdaStateDependent(const VectorXd& new_lambda) {

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
      stdMatrixXd lambdas = TensorPhyloUtils::VectorToStdMat(new_lambda, 1);

      // set the value
      internal->setLambda(lambda_times, lambdas);

    }

    // set time-dependent lambda
    void setLambdaTimeDependent(const VectorXd& new_lambda_times, const VectorXd& new_lambda) {

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

      // copy the time-dependent rates per state
      stdMatrixXd lambdas(new_lambda.size());
      for(size_t i = 0; i < new_lambda.size(); ++i) {
        lambdas.at(i) = stdVectorXd(dim, new_lambda(i));
      }

      // set the value
      internal->setLambda(lambda_times, lambdas);

    }

    // set time/state dependent lambda
    void setLambdaTimeStateDependent(const VectorXd& new_lambda_times, const MatrixXd& new_lambda) {

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

    void setMuConstant(const double& new_mu) {

      if ( safe ) {

        // check that it's valid rate
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_mu) == false ) {
          stop("Error setting extinction rates. Rates must be strictly non-negative.");
        }

      }

      // set mu times to empty
      stdVectorXd mu_times = stdVectorXd();

      // repeat the mus
      stdMatrixXd mus = TensorPhyloUtils::ScalarToStdMat(new_mu, dim);

      // set the value
      internal->setMu(mu_times, mus);

    }

    // set state dependent
    void setMuStateDependent(const VectorXd& new_mu) {

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
      stdMatrixXd mus = TensorPhyloUtils::VectorToStdMat(new_mu, 1);

      // set the value
      internal->setMu(mu_times, mus);

    }

    // set time-dependent mu
    void setMuTimeDependent(const VectorXd& new_mu_times, const VectorXd& new_mu) {

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

      // copy the time-dependent rates per state
      stdMatrixXd mus(new_mu.size());
      for(size_t i = 0; i < new_mu.size(); ++i) {
        mus.at(i) = stdVectorXd(dim, new_mu(i));
      }

      // set the value
      internal->setMu(mu_times, mus);

    }

    // set time/state dependent mu
    void setMuTimeStateDependent(const VectorXd& new_mu_times, const MatrixXd& new_mu) {

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

    void setPhiConstant(const double& new_phi) {

      if ( safe ) {

        // check that it's valid rate
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
          stop("Error setting sampling rates. Rates must be strictly non-negative.");
        }

      }

      // set times to empty
      stdVectorXd phi_times = stdVectorXd();

      // repeat the phis
      stdMatrixXd phis = TensorPhyloUtils::ScalarToStdMat(new_phi, dim);

      // set the value
      internal->setPhi(phi_times, phis);

    }

    // set state dependent
    void setPhiStateDependent(const VectorXd& new_phi) {

      if ( safe ) {

        // check that rates are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_phi) == false ) {
          stop("Error setting sampling rates. Rates must be strictly non-negative.");
        }

        // check the size
        if ( TensorPhyloUtils::hasDimensions(new_phi, dim) == false ) {
          stop("Error setting sampling rates. Number of rates must equal the number of states.");
        }

      }

      // set times to empty
      stdVectorXd phi_times = stdVectorXd();

      // create the rates
      stdMatrixXd phis = TensorPhyloUtils::VectorToStdMat(new_phi, 1);

      // set the value
      internal->setPhi(phi_times, phis);

    }

    // set time-dependent phi
    void setPhiTimeDependent(const VectorXd& new_phi_times, const VectorXd& new_phi) {

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

      // copy the time-dependent rates per state
      stdMatrixXd phis(new_phi.size());
      for(size_t i = 0; i < new_phi.size(); ++i) {
        phis.at(i) = stdVectorXd(dim, new_phi(i));
      }

      // set the value
      internal->setPhi(phi_times, phis);

    }

    // set time/state dependent phi
    void setPhiTimeStateDependent(const VectorXd& new_phi_times, const MatrixXd& new_phi) {

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

    void setDeltaConstant(const double& new_delta) {

      if ( safe ) {

        // check that it's valid rate
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
          stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
        }

      }

      // set delta times to empty
      stdVectorXd delta_times = stdVectorXd();

      // repeat the deltas
      stdMatrixXd deltas = TensorPhyloUtils::ScalarToStdMat(new_delta, dim);

      // set the value
      internal->setDelta(delta_times, deltas);

    }

    // set state dependent
    void setDeltaStateDependent(const VectorXd& new_delta) {

      if ( safe ) {

        // check that it's valid rate
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_delta) == false ) {
          stop("Error setting destructive-sampling rates. Rates must be strictly non-negative.");
        }

        // check the size
        if ( TensorPhyloUtils::hasDimensions(new_delta, dim) == false ) {
          stop("Error setting destructive-sampling rates. Number of rates must equal the number of states.");
        }

      }

      // set delta times to empty
      stdVectorXd delta_times = stdVectorXd();

      // create the rates
      stdMatrixXd deltas = TensorPhyloUtils::VectorToStdMat(new_delta, 1);

      // set the value
      internal->setDelta(delta_times, deltas);

    }

    // set time-dependent delta
    void setDeltaTimeDependent(const VectorXd& new_delta_times, const VectorXd& new_delta) {

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

      // copy the time-dependent rates per state
      stdMatrixXd deltas(new_delta.size());
      for(size_t i = 0; i < new_delta.size(); ++i) {
        deltas.at(i) = stdVectorXd(dim, new_delta(i));
      }

      // set the value
      internal->setDelta(delta_times, deltas);

    }

    // set time/state dependent delta
    void setDeltaTimeStateDependent(const VectorXd& new_delta_times, const MatrixXd& new_delta) {

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
        if ( TensorPhyloUtils::hasDimensions(new_delta, new_delta_times.size() + 1, dim) == false ) {
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

    void setEtaConstantEqual(const double& new_eta) {

      if ( safe ) {

        // check that it's valid rate
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta) == false ) {
          stop("Error setting transition rates. Rates must be strictly non-negative.");
        }

      }

      // set delta times to empty
      stdVectorXd eta_times = stdVectorXd();

      // create a matrix
      MatrixXd tmp       = MatrixXd::Constant(dim, dim, new_eta) / ((double)dim - 1);
      tmp.diagonal()     = VectorXd::Constant(dim, -new_eta);
      stdMatrixXd tmpStd = TensorPhyloUtils::EigenToStd(tmp);

      // create the vector of etas
      std::vector<stdMatrixXd> etas = std::vector<stdMatrixXd>(1, tmpStd);

      // set the value
      internal->setEta(eta_times, etas);

    }

    void setEtaConstantUnequal(const RateMatrix& new_eta) {

      if ( safe ) {

        // make sure this is a rate matrix
        if ( TensorPhyloUtils::isRateMatrix(new_eta) == false ) {
          stop("Error setting transition rates. Rate matrix must be square, and the rows must sum to zero.");
        }

        // check the dimensions
        if ( TensorPhyloUtils::hasDimensions(new_eta, dim, dim) == false ) {
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

    void setEtaTimeDependentEqual(const VectorXd& new_eta_times, const VectorXd& new_eta) {

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
        MatrixXd tmp   = MatrixXd::Constant(dim, dim, new_eta[i]) / ((double)dim - 1);
        tmp.diagonal() = VectorXd::Constant(dim, -new_eta[i]);
        etas.at(i) = TensorPhyloUtils::EigenToStd(tmp);
      }

      // set the value
      internal->setEta(eta_times, etas);

    }

    void setEtaTimeDependentUnequal(const VectorXd& new_eta_times, const RateMatrixList& new_eta) {

      if ( safe ) {

        // check that times are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_eta_times) == false ) {
          stop("Error setting transition rates. Times must be strictly non-negative.");
        }

        // check that there are the right number of matrices
        if ( (new_eta_times.size() + 1) != new_eta.size() ) {
          stop("Error setting transition rates. Number of change times must be 1 less than the number of rate matrices.");
        }

        // loop over each matrix
        for(size_t i = 0; i < new_eta.size(); ++i) {

          // make sure this is a rate matrix
          if ( TensorPhyloUtils::isRateMatrix( new_eta.at(i) ) == false ) {
            stop("Error setting transition rates. Rate matrix must be square, and the rows must sum to zero.");
          }

          // check the dimensions
          if ( TensorPhyloUtils::hasDimensions(new_eta.at(i), dim, dim) == false ) {
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

    // KNOWN BUG: if you call computeLogLikelihood on a model without omega, then set omega,
    // the scheduler does not correctly recognize that it needs to update.
    void setOmegaConstant(const CladoEvents& new_omega) {

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

      // workaround for issue with omega
      // force an update to the scheduler
      if ( has_omega == false ) {
        internal->forceSchedulerUpdate();
        internal->forceApproximatorDirty();
        has_omega = true;
      }

    }

    void setOmegaTimeDependent(const VectorXd& new_omega_times, const CladoEventsList& new_omegas) {

      if ( safe ) {

        // check that times are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_omega_times) == false ) {
          stop("Error setting cladogenetic events. Times must be strictly non-negative.");
        }

        // check that there are the right number of cladogenetic arrays
        if ( (new_omega_times.size() + 1) != new_omegas.size() ) {
          stop("Error setting cladogenetic events. Number of change times must be 1 less than the number of cladogenetic arrays.");
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

      // workaround for issue with omega
      // force an update to the scheduler
      if ( has_omega == false ) {
        internal->forceSchedulerUpdate();
        internal->forceApproximatorDirty();
        has_omega = true;
      }

    }

    /////////////////////
    // mass speciation //
    /////////////////////

    void setUpsilonConstant(const VectorXd& new_upsilon_times, const VectorXd& new_upsilon) {

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

      // copy the time-dependent rates per state
      stdMatrixXd upsilons(new_upsilon.size());
      for(size_t i = 0; i < new_upsilon.size(); ++i) {
        upsilons.at(i) = stdVectorXd(dim, new_upsilon(i));
      }

      // set value
      internal->setMassSpeciationEvents(upsilon_times, upsilons);

    }

    void setUpsilonStateDependent(const VectorXd& new_upsilon_times, const MatrixXd& new_upsilon) {

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

    void setGammaConstant(const VectorXd& new_gamma_times, const VectorXd& new_gamma) {

      if ( safe ) {

        // check that times are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
          stop("Error setting mass-extinction events. Times must be strictly non-negative.");
        }

        // make sure the number of times and magnitudes match
        if ( TensorPhyloUtils::hasDimensions(new_gamma, new_gamma_times.size()) == false ) {
          stop("Error setting mass-extinction events. Number of times must equal to the number of magnitudes.");
        }

        // make sure the gammas are probabilities
        if ( TensorPhyloUtils::isProbability(new_gamma) == false ) {
          stop("Error setting mass-extinction events. Magnitudes must be between 0 and 1 (inclusive).");
        }

      }

      // set the time variable
      stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

      // copy the time-dependent rates per state
      stdMatrixXd gammas(new_gamma.size());
      for(size_t i = 0; i < new_gamma.size(); ++i) {
        gammas.at(i) = stdVectorXd(dim, new_gamma(i));
      }

      // set zeta to identity
      MatrixXd tmp( MatrixXd::Identity(dim, dim) );
      stdMatrixXd zeta = TensorPhyloUtils::EigenToStd( tmp );
      std::vector<stdMatrixXd> zetas( new_gamma_times.size(), zeta);

      // set value
      internal->setMassExtinctionEvents(gamma_times, gammas);
      internal->setMassExtinctionStateChangeProb(zetas);

    }

    void setGammaStateDependent(const VectorXd& new_gamma_times, const MatrixXd& new_gamma) {

      if ( safe ) {

        // check that times are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
          stop("Error setting mass-extinction events. Times must be strictly non-negative.");
        }

        // make sure the number of times and magnitudes match
        if ( TensorPhyloUtils::hasDimensions(new_gamma, new_gamma_times.size(), dim) == false ) {
          stop("Error setting mass-speciation events. Number of times must equal to the number of magnitudes.");
        }

        // make sure the gammas are probabilities
        if ( TensorPhyloUtils::isProbability(new_gamma) == false ) {
          stop("Error setting mass-speciation events. Magnitudes must be between 0 and 1 (inclusive).");
        }

      }

      // set the time variable
      stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

      // set the rate variable
      stdMatrixXd gammas = TensorPhyloUtils::EigenToStd(new_gamma);

      // set zeta to identity
      MatrixXd tmp( MatrixXd::Identity(dim, dim) );
      stdMatrixXd zeta = TensorPhyloUtils::EigenToStd( tmp );
      std::vector<stdMatrixXd> zetas( new_gamma_times.size(), zeta);

      // set value
      internal->setMassExtinctionEvents(gamma_times, gammas);
      internal->setMassExtinctionStateChangeProb(zetas);

    }

    void setGammaAndZetaConstant(const VectorXd& new_gamma_times, const VectorXd& new_gamma, const ProbabilityMatrixList& new_zeta) {

      if ( safe ) {

        // check that times are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
          stop("Error setting mass-extinction events. Times must be strictly non-negative.");
        }

        // make sure the number of times and magnitudes match
        if ( TensorPhyloUtils::hasDimensions(new_gamma, new_gamma_times.size()) == false ) {
          stop("Error setting mass-extinction events. Number of times must equal to the number of magnitudes.");
        }

        // make sure the gammas are probabilities
        if ( TensorPhyloUtils::isProbability(new_gamma) == false ) {
          stop("Error setting mass-extinction events. Magnitudes must be between 0 and 1 (inclusive).");
        }

        // check zeta
        if ( new_gamma_times.size() != new_zeta.size() ) {
          stop("Error setting mass-extinction state-change events. Number of times must equal to the number of transition probability matrices.");
        }

        // loop over each element of list
        for(size_t i = 0; i < new_zeta.size(); ++i) {

          // check dimensions of matrix
          if ( TensorPhyloUtils::hasDimensions(new_zeta.at(i), dim, dim) == false ) {
            stop("Error setting mass-extinction state-change events. Transition probability matrix must be square.");
          }

          // check that it's a transition probability matrix
          if ( TensorPhyloUtils::isTransitionProbabilityMatrix(new_zeta.at(i)) == false ) {
            stop("Error setting mass-extinction state-change events. Rows of transition probability matrix must sum to 1.0.");
          }

        }

      }

      // set stuff
      stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

      // copy the time-dependent rates per state
      stdMatrixXd gammas(new_gamma.size());
      for(size_t i = 0; i < new_gamma.size(); ++i) {
        gammas.at(i) = stdVectorXd(dim, new_gamma(i));
      }

      // copy the transition-probability matrices
      std::vector<stdMatrixXd> zetas( new_zeta.size() );
      for(size_t i = 0; i < new_zeta.size(); ++i) {
        zetas.at(i) = TensorPhyloUtils::EigenToStd(new_zeta.at(i).getMatrix());
      }

      // set value
      internal->setMassExtinctionEvents(gamma_times, gammas);
      internal->setMassExtinctionStateChangeProb(zetas);

    }

    void setGammaAndZetaStateDependent(const VectorXd& new_gamma_times, const MatrixXd& new_gamma, const ProbabilityMatrixList& new_zeta) {

        if ( safe ) {

          // check that times are valid
          if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
            stop("Error setting mass-extinction events. Times must be strictly non-negative.");
          }

          // make sure the number of times and magnitudes match
          if ( TensorPhyloUtils::hasDimensions(new_gamma, new_gamma_times.size(), dim) == false ) {
            stop("Error setting mass-speciation events. Number of times must equal to the number of magnitudes.");
          }

          // make sure the gammas are probabilities
          if ( TensorPhyloUtils::isProbability(new_gamma) == false ) {
            stop("Error setting mass-speciation events. Magnitudes must be between 0 and 1 (inclusive).");
          }

          // check zeta
          if ( new_gamma_times.size() != new_zeta.size() ) {
            stop("Error setting mass-extinction state-change events. Number of times must equal to the number of transition probability matrices.");
          }

          // loop over each element of list
          for(size_t i = 0; i < new_zeta.size(); ++i) {

            // check dimensions of matrix
            if ( TensorPhyloUtils::hasDimensions(new_zeta.at(i), dim, dim) == false ) {
              stop("Error setting mass-extinction state-change events. Transition probability matrix must be square.");
            }

            // check that it's a transition probability matrix
            if ( TensorPhyloUtils::isTransitionProbabilityMatrix(new_zeta.at(i)) == false ) {
              stop("Error setting mass-extinction state-change events. Rows of transition probability matrix must sum to 1.0.");
            }

          }

        }

        // set the time variable
        stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

        // set the rate variable
        stdMatrixXd gammas = TensorPhyloUtils::EigenToStd(new_gamma);

        // copy the transition-probability matrices
        std::vector<stdMatrixXd> zetas( new_zeta.size() );
        for(size_t i = 0; i < new_zeta.size(); ++i) {
          zetas.at(i) = TensorPhyloUtils::EigenToStd(new_zeta.at(i).getMatrix());
        }

        // set value
        internal->setMassExtinctionEvents(gamma_times, gammas);
        internal->setMassExtinctionStateChangeProb(zetas);

    }

    /////////////////////////////////////////////////
    // mass state-change (without mass-extinction) //
    /////////////////////////////////////////////////

    void setZeta(const VectorXd& new_gamma_times, const ProbabilityMatrixList& new_zeta) {

      if ( safe ) {

        // check that times are valid
        if ( TensorPhyloUtils::isStrictlyNonNegative(new_gamma_times) == false ) {
          stop("Error setting mass-state-change events. Times must be strictly non-negative.");
        }

        // check zeta
        if ( new_gamma_times.size() != new_zeta.size() ) {
          stop("Error setting mass-state-change events. Number of times must equal to the number of transition probability matrices.");
        }

        // loop over each element of list
        for(size_t i = 0; i < new_zeta.size(); ++i) {

          // check dimensions of matrix
          if ( TensorPhyloUtils::hasDimensions(new_zeta.at(i), dim, dim) == false ) {
            stop("Error setting mass-state-change events. Transition probability matrix must be square.");
          }

          // check that it's a transition probability matrix
          if ( TensorPhyloUtils::isTransitionProbabilityMatrix(new_zeta.at(i)) == false ) {
            stop("Error setting mass-state-change events. Rows of transition probability matrix must sum to 1.0.");
          }

        }

      }

      // set stuff
      stdVectorXd gamma_times = TensorPhyloUtils::EigenToStd(new_gamma_times);

      // create some empty mass-extinction events
      stdMatrixXd gammas( new_gamma_times.size(), stdVectorXd(dim, 0));

      // copy the transition-probability matrices
      std::vector<stdMatrixXd> zetas( new_zeta.size() );
      for(size_t i = 0; i < new_zeta.size(); ++i) {
        zetas.at(i) = TensorPhyloUtils::EigenToStd(new_zeta.at(i).getMatrix());
      }

      // set value
      internal->setMassExtinctionEvents(gamma_times, gammas);
      internal->setMassExtinctionStateChangeProb(zetas);

    }

    ///////////////////
    // mass sampling //
    ///////////////////

    void setRhoPresent(const double& new_rho) {

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

    void setRhoPresentStateDependent(const VectorXd& new_rho) {

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

    void setRhoConstant(const VectorXd& new_rho_times, const VectorXd& new_rho) {

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

      // copy the time-dependent rates per state
      stdMatrixXd rhos(new_rho.size());
      for(size_t i = 0; i < new_rho.size(); ++i) {
        rhos.at(i) = stdVectorXd(dim, new_rho(i));
      }

      // set value
      internal->setMassSamplingEvents(rho_times, rhos);

    }

    void setRhoStateDependent(const VectorXd& new_rho_times, const MatrixXd& new_rho) {

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

    void setXiConstant(const VectorXd& new_xi_times, const VectorXd& new_xi) {

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

      // copy the time-dependent rates per state
      stdMatrixXd xis(new_xi.size());
      for(size_t i = 0; i < new_xi.size(); ++i) {
        xis.at(i) = stdVectorXd(dim, new_xi(i));
      }

      // set value
      internal->setMassDestrSamplingEvents(xi_times, xis);

    }

    void setXiStateDependent(const VectorXd& new_xi_times, const MatrixXd& new_xi) {

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

    ///////////////////////////////////
    // the quasistationary frequency //
    ///////////////////////////////////

    VectorXd getQuasiStationaryFrequency(double t) {

      // get the frequency
      VectorXd qsf = internal->getQuasiStationaryFrequency(t);

      // return it as a vectorXd
      return qsf;

    }

  private:

    // init
    void init() {

      if ( ready ) {
        return;
      }

      // create the internal
      internal = DistributionHandlerImpl::create();

      // some default things
      internal->setConditionalProbCompatibilityMode(false);
    	internal->setNumberOfThreads(1);
      internal->setConditionalProbabilityType(conditionalProbability_t::TIME);
      internal->setIntegrationScheme(integrationScheme_t::RUNGE_KUTTA54);

      // default root frequency is flat
      setRootPriorFlat();

      // default speciation rate is 1.0
      setLambdaConstant(1.0);

      // default extinction rate is 0.0
      setMuConstant(0.0);

      // default sampling rates are zero
      setPhiConstant(0.0);
      setDeltaConstant(0.0);

      // default transition rate is 1.0
      setEtaConstantEqual(1.0);

      // default mass-sampling probability is 1 at the present
      setRhoPresent(1.0);

      // parameters not set/not in default model:
      // omega
      // upsilon
      // gamma nor zeta
      // xi

    }

    // pointer
    boost::shared_ptr<DistributionHandlerImpl> internal;

    // safe mode
    bool ready;
    bool safe;

    // dimensions
    size_t dim;
    std::vector<std::string> state_labels;

    // the tree
    List phy;

    // workaround for omega problem
    bool has_omega;

};

#endif
