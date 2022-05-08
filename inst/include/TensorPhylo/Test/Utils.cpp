/*
 * Utils.cpp
 *
 *  Created on: Sep 16, 2019
 *      Author: xaviermeyer
 */

#include "Utils.h"

#include <iostream>
#include <map>

#include "../Likelihood/Approximator/AsynchronousApproximator.h"
#include "../Likelihood/Approximator/BaseApproximator.h"
#include "../Likelihood/Approximator/Specialized/DenseUnobserved.h"
#include "../Likelihood/Approximator/StochasticMapping/SequentialStochasticMappingCPU.h"
#include "../Likelihood/Kernels/CPU/EigenUtilsMatrix.h"
#include "../Likelihood/Scheduler/DAG/DAG.h"
#include "Interface/DistributionHandler.h"
#include "Interface/DistributionHandlerImpl.h"
#include "Utils/Profiling/CustomProfiling.h"
#include "Utils/Profiling/UniqueProfiler.h"
#include "Likelihood/Scheduler/DAG/IncDAG.h"

namespace Test {

/****************************************************************************************************************/
/**********************************************    DEFINITIONS    ***********************************************/
/****************************************************************************************************************/

namespace Definitions {
	//static const std::vector< std::string > intSchemeNames = boost::assign::list_of("EULER")("RUNGE_KUTTA4")("RUNGE_KUTTA_DOPRI5")("RUNGE_KUTTA_DOPRI5");
}

/****************************************************************************************************************/
/****************************************    INIT DATA FUNCTIONS    ********************************************/
/****************************************************************************************************************/

namespace InitData {

Phylogeny::NexusReader::NexusParserSharedPtr nexusFromSingleBranchFile(std::string fileName) {

	std::ifstream iFile(fileName.c_str(), std::ios::in);
	assert(iFile.good() && "Single branch input file not found.");

	size_t nCat = 0;
	size_t nNode = 0;

	iFile >> nCat;
	iFile >> nNode;

	std::vector<double> times;
	std::vector< std::vector<double> > initStates;

	for(size_t iN=0; iN<nNode; ++iN) {
		double time = 0.;
		iFile >> time;
		times.push_back(time);

		if(iN == nNode -1) break;

		std::vector<double> initState(nCat);
		for(size_t iS=0; iS<nCat; ++iS) {
			iFile >> initState[iS];
		}
		initStates.push_back(initState);
	}

	Phylogeny::NexusReader::NexusParserSharedPtr nexusParser( new Phylogeny::NexusReader::NexusParser(nCat, times, initStates) );

	return nexusParser;
}

Phylogeny::NexusReader::NexusParserSharedPtr nexusFromFile(std::string fileName) {
	Phylogeny::NexusReader::NexusParserSharedPtr nexusParser( new Phylogeny::NexusReader::NexusParser(fileName) );
	return nexusParser;
}

Phylogeny::NexusReader::NexusParserSharedPtr nexusFromFile(std::string nexusFilename, std::string tsvFilename) {
	Phylogeny::NexusReader::NexusParserSharedPtr nexusParser( new Phylogeny::NexusReader::NexusParser(nexusFilename, tsvFilename) );
	return nexusParser;
}

Parameters::ContainerSharedPtr parametersFromFile(std::string fileName) {
	Parameters::ContainerSharedPtr parameters(new Parameters::Container);
	parameters->readFromFile(fileName);
	return parameters;
}

Parameters::ContainerSharedPtr parametersYule1(size_t nState, double lambda, double eta, double rho) {
	Parameters::ContainerSharedPtr parameters( new Parameters::Container );

	size_t K = nState;

	parameters->intLikApproximator = 1;
	parameters->applyTreeLikCorrection = false;

	parameters->intScheme = 3;
	parameters->condType = 0;

	parameters->deltaT = 0.001;

	parameters->rootPrior = (1./K)*Eigen::VectorXd::Ones(K);

	parameters->lambda = lambda*Eigen::VectorXd::Ones(K);

	parameters->mu = Eigen::VectorXd::Ones(K);
	parameters->mu(0) = 0.0;
	parameters->mu(1) = 0.0;

	parameters->phi   = Eigen::VectorXd::Zero(K);
	parameters->delta = Eigen::VectorXd::Zero(K);

	parameters->massSamplingTimes.push_back(0.0);
	Eigen::VectorXd RhoAtZero = rho * Eigen::VectorXd::Ones(K);
	parameters->massSamplingProb.push_back(RhoAtZero);

	parameters->eta.resize(K,K);
	for(size_t i=0; i<K; ++i) {
		for(size_t j=0; j<K; ++j) {
				if(i == j) parameters->eta(i,j) = -(K-1.0) * eta;
				else parameters->eta(i,j) = eta;
		}
	}

	parameters->omega.resize(nState);
	for(size_t iX=0; iX<nState; ++iX) {
		parameters->omega[iX].resize(nState, nState);
		parameters->omega[iX].setZero();
		for(size_t iY=0; iY<nState; ++iY) {
			for(size_t iZ=0; iZ<nState; ++iZ) {
				if(iX == iY && iY == iZ) {
					parameters->omega[iX](iX,iY) = 1.0;
				}
			}
		}
	}

	return parameters;
}

Parameters::ContainerSharedPtr parametersBD1(size_t nState, double lambda, double mu, double eta, double rho) {
	Parameters::ContainerSharedPtr parameters( new Parameters::Container );

	size_t K = nState;

	parameters->intLikApproximator = 1;
	parameters->applyTreeLikCorrection = false;

	parameters->intScheme = 3;
	parameters->condType = 0;

	parameters->deltaT = 0.001;

	parameters->rootPrior = (1.0/K) * Eigen::VectorXd::Ones(K);

	parameters->lambda = lambda * Eigen::VectorXd::Ones(K);

	parameters->mu = mu * Eigen::VectorXd::Ones(K);

	parameters->phi   = Eigen::VectorXd::Zero(K);
	parameters->delta = Eigen::VectorXd::Zero(K);

	parameters->massSamplingTimes.push_back(0.0);
	Eigen::VectorXd RhoAtZero = rho * Eigen::VectorXd::Ones(K);
	parameters->massSamplingProb.push_back(RhoAtZero);

	parameters->eta.resize(K,K);
	for(size_t i=0; i<K; ++i) {
		for(size_t j=0; j<K; ++j) {
				if(i == j) parameters->eta(i,j) = -(K-1.0) * eta;
				else parameters->eta(i,j) = eta;
		}
	}

	parameters->omega.resize(nState);
	for(size_t iX=0; iX<nState; ++iX) {
		parameters->omega[iX].resize(nState, nState);
		parameters->omega[iX].setZero();
		for(size_t iY=0; iY<nState; ++iY) {
			for(size_t iZ=0; iZ<nState; ++iZ) {
				if(iX == iY && iY == iZ) {
					parameters->omega[iX](iX,iY) = 1.0;
				}
			}
		}
	}

	return parameters;
}

Parameters::ContainerSharedPtr parametersFBD1(size_t nState, double lambda, double mu, double phi, double eta, double rho) {
	Parameters::ContainerSharedPtr parameters( new Parameters::Container );

	size_t K = nState;

	parameters->intLikApproximator = 1;
	parameters->applyTreeLikCorrection = false;

	parameters->intScheme = 3;
	parameters->condType = 0;

	parameters->deltaT = 0.001;

	parameters->rootPrior = (1.0/K) * Eigen::VectorXd::Ones(K);

	parameters->lambda = lambda * Eigen::VectorXd::Ones(K);
	parameters->mu     = mu * Eigen::VectorXd::Ones(K);
	parameters->phi    = phi * Eigen::VectorXd::Ones(K);
	parameters->delta = Eigen::VectorXd::Zero(K);

	parameters->massSamplingTimes.push_back(0.0);
	Eigen::VectorXd RhoAtZero = rho * Eigen::VectorXd::Ones(K);
	parameters->massSamplingProb.push_back(RhoAtZero);

	parameters->eta.resize(K,K);
	for(size_t i=0; i<K; ++i) {
		for(size_t j=0; j<K; ++j) {
				if(i == j) parameters->eta(i,j) = -(K-1.0) * eta;
				else parameters->eta(i,j) = eta;
		}
	}

	parameters->omega.resize(nState);
	for(size_t iX=0; iX<nState; ++iX) {
		parameters->omega[iX].resize(nState, nState);
		parameters->omega[iX].setZero();
		for(size_t iY=0; iY<nState; ++iY) {
			for(size_t iZ=0; iZ<nState; ++iZ) {
				if(iX == iY && iY == iZ) {
					parameters->omega[iX](iX,iY) = 1.0;
				}
			}
		}
	}

	return parameters;
}


Parameters::ContainerSharedPtr parametersBDME(size_t nState, double lambda, double mu, double eta, double rho, double me_time, double me_prob) {
	Parameters::ContainerSharedPtr parameters( new Parameters::Container );

	size_t K = nState;

	parameters->intLikApproximator = 1;
	parameters->applyTreeLikCorrection = false;

	parameters->intScheme = 3;
	parameters->condType = 0;

	parameters->deltaT = 0.001;

	parameters->rootPrior = (1.0/K) * Eigen::VectorXd::Ones(K);

	parameters->lambda = lambda * Eigen::VectorXd::Ones(K);

	parameters->mu = mu * Eigen::VectorXd::Ones(K);

	parameters->phi   = Eigen::VectorXd::Zero(K);
	parameters->delta = Eigen::VectorXd::Zero(K);

	parameters->massSamplingTimes.push_back(0.0);
	Eigen::VectorXd RhoAtZero = rho * Eigen::VectorXd::Ones(K);
	parameters->massSamplingProb.push_back(RhoAtZero);

	parameters->massExtinctionTimes.push_back(me_time);
	Eigen::VectorXd GammaAtZero = me_prob * Eigen::VectorXd::Ones(K);
	parameters->massExtinctionProb.push_back(GammaAtZero);
	Eigen::MatrixXd Zeta = Eigen::MatrixXd::Identity(K, K);
	parameters->massExtinctionStateChangeProb.push_back(Zeta);

	parameters->eta.resize(K,K);
	for(size_t i=0; i<K; ++i) {
		for(size_t j=0; j<K; ++j) {
				if(i == j) parameters->eta(i,j) = -(K-1.0) * eta;
				else parameters->eta(i,j) = eta;
		}
	}

	parameters->omega.resize(nState);
	for(size_t iX=0; iX<nState; ++iX) {
		parameters->omega[iX].resize(nState, nState);
		parameters->omega[iX].setZero();
		for(size_t iY=0; iY<nState; ++iY) {
			for(size_t iZ=0; iZ<nState; ++iZ) {
				if(iX == iY && iY == iZ) {
					parameters->omega[iX](iX,iY) = 1.0;
				}
			}
		}
	}

	return parameters;
}

Parameters::ContainerSharedPtr parametersBiSSE(size_t nState, std::vector<double> lambda, std::vector<double> mu, double eta, double rho) {

	Parameters::ContainerSharedPtr parameters( new Parameters::Container );

	size_t K = nState;

	parameters->intLikApproximator = 1;
	parameters->applyTreeLikCorrection = false;

	parameters->intScheme = 3;
	parameters->condType = 0;

	parameters->deltaT = 0.001;

	parameters->rootPrior = (1.0/K) * Eigen::VectorXd::Ones(K);

	parameters->lambda = Eigen::VectorXd::Zero(K);
	parameters->mu     = Eigen::VectorXd::Zero(K);
	for(size_t i = 0; i < K; ++i) {
		parameters->lambda(i) = lambda[i];
		parameters->mu(i)     = mu[i];
	}

	parameters->phi   = Eigen::VectorXd::Zero(K);
	parameters->delta = Eigen::VectorXd::Zero(K);

	parameters->massSamplingTimes.push_back(0.0);
	Eigen::VectorXd RhoAtZero = rho * Eigen::VectorXd::Ones(K);
	parameters->massSamplingProb.push_back(RhoAtZero);

	parameters->eta.resize(K,K);
	for(size_t i=0; i<K; ++i) {
		for(size_t j=0; j<K; ++j) {
				if(i == j) parameters->eta(i,j) = -(K-1.0) * eta;
				else parameters->eta(i,j) = eta;
		}
	}

	parameters->omega.resize(nState);
	for(size_t iX=0; iX<nState; ++iX) {
		parameters->omega[iX].resize(nState, nState);
		parameters->omega[iX].setZero();
		for(size_t iY=0; iY<nState; ++iY) {
			for(size_t iZ=0; iZ<nState; ++iZ) {
				if(iX == iY && iY == iZ) {
					parameters->omega[iX](iX,iY) = 1.0;
				}
			}
		}
	}

	return parameters;

}


} /* namespace InitData */

/****************************************************************************************************************/
/************************************    SEQUENTIAL CPU FUNCTIONS    ***************************************/
/****************************************************************************************************************/
namespace SequentialCPU {

double computeLikForTests(Likelihood::Integrator::integrationScheme_t intScheme,
				  Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
				  Parameters::ContainerSharedPtr parameters) {

	double lLik = 0.;

	// Making sure we test parallel version in parallel (2 therads min)
#ifdef _OPENMP
	if((parameters->intLikApproximator == TensorPhylo::Interface::PARALLEL_BRANCHWISE ||
		parameters->intLikApproximator == TensorPhylo::Interface::PARALLEL_OPTIMIZED) &&
			parameters->nThreads == 1) {

		parameters->nThreads = 2;

		lLik = computeLogLik(1, intScheme, nexusParser, parameters);

		parameters->nThreads = 1;

	} else {
		lLik = computeLogLik(1, intScheme, nexusParser, parameters);
	}
#else
	lLik = computeLogLik(1, intScheme, nexusParser, parameters);
#endif

	return lLik;
}

double computeLikForInterface(Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
							  const std::string &parametersFile) {

	//Utils::Parallel::Manager::getInstance()->setNThread(1);
	//Utils::Parallel::Manager::getInstance()->setMaxNThread(1);

	boost::shared_ptr<TensorPhylo::Interface::DistributionHandlerImpl> ptrDH(TensorPhylo::Interface::DistributionHandlerImpl::create());
	//ptrDH->setDebugMode(TensorPhylo::Interface::DBG_PRINT);

	// *********** WARNING ORDER IS IMPORTANT : LoadStateFromFile MUST BE CALLED FIRST! **************/

	// Loading the parameters
	ptrDH->loadStateFromFile(parametersFile);

	// Loading tree:
	std::vector<std::string> taxaNames = nexusParser->getTaxaNames();
	std::map<int, std::string> taxaIdToNameMap(nexusParser->getTaxaIdToNameMap());
	PS::TreeSharedPtr ptrTree( new PS::Tree(nexusParser) );
	std::string newickStr(ptrTree->getNewickString(taxaIdToNameMap));
	ptrDH->setTree(newickStr);

	// Trick to get probes automatically spaced
	if(ptrDH->synchMonitoring.size() == 1 && ptrDH->synchMonitoring.front() < 0.) {
		double step = fabs(ptrDH->synchMonitoring.front());
		ptrDH->synchMonitoring.clear();
		double ageOfOldest = ptrTree->getOldestNode()->getAge();

		for(double time=1.e-5; time < ageOfOldest; time += step) {
			ptrDH->synchMonitoring.push_back(time);
		}
		ptrDH->synchMonitoring.push_back(ageOfOldest-1.e-5);
	}

	// Loading data:
	if(nexusParser->hasData()) {
		Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

		// Loading the data
		std::map<std::string, std::vector<double> > probabilityMap;
		for(size_t iT=0; iT<taxaNames.size(); ++iT) {
			std::vector<double> prob(nexusParser->getNState(), 0.);
			// Load
			Eigen::VectorXd pVec = ptrData->getProbForTaxaLabel(taxaNames[iT]);
			// To std vector
			for(size_t iP=0; iP<prob.size(); ++iP) prob[iP] = pVec(iP);
			// in map
			probabilityMap[taxaNames[iT]] = prob;
		}

		ptrDH->setData(nexusParser->getTaxaNames(), probabilityMap);
	}

	Utils::Parallel::Manager::getInstance()->setMaxNThread(ptrDH->nThreads);
	Utils::Parallel::Manager::getInstance()->setNThread(ptrDH->nThreads);

	double logLik = ptrDH->computeLogLikelihood();

	if(Utils::Output::outputManager().check(Utils::Output::LOW_VERB) || Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB) ) {
		std::cout << "LogLik = " << logLik << std::scientific << std::setprecision(8) << " ( scientific: " << logLik << " ) " <<  std::endl;
	}

	if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB)) {
		ptrDH->approximator->orderProbes();
		const std::vector<Likelihood::Monitor::ProbeState>& probes = ptrDH->approximator->getObservedProbesState();
		std::cout << "------- PROBES --------" << std::endl;
		for(size_t i=0; i<probes.size(); ++i) {
			std::cout << probes[i].toString() << std::endl;
		}
	}

	if(Utils::Output::outputManager().check(Utils::Output::DUMP_VERB)) {
		const std::vector<Likelihood::Monitor::ProbeState>& probes = ptrDH->approximator->getObservedProbesState();
		std::cout << "------- DUMPING PROBES IN 'dump_probes.txt' --------" << std::endl;
		std::ofstream oFile(Utils::Output::outputManager().getDumpFileName());
		for(size_t i=0; i<probes.size(); ++i) {
			oFile << probes[i].toDumpFormat();
		}

		// Here we recompute the likelihood without probes but we keep track of interation times (t+dt)
		std::vector<double> emptyVec;
		ptrDH->setSyncMonitors(emptyVec);
		logLik = ptrDH->computeLogLikelihood();
		const std::vector<double> &intTimes = ptrDH->approximator->getIntegrationTimes();
		oFile << "Integration times" << std::endl;
		for(size_t i=0; i<intTimes.size(); ++i) {
			oFile << intTimes[i] << " ";
		}
		oFile << std::endl;
	}

	return logLik;

}

TensorPhylo::Interface::vecHistories_t computeStochasticMappingForInterface(size_t nHistories,
											Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
							  	  	  	  	const std::string &parametersFile) {

	//Utils::Parallel::Manager::getInstance()->setNThread(1);
	//Utils::Parallel::Manager::getInstance()->setMaxNThread(1);

	boost::shared_ptr<TensorPhylo::Interface::DistributionHandlerImpl> ptrDH(TensorPhylo::Interface::DistributionHandlerImpl::create());
	//ptrDH->setDebugMode(TensorPhylo::Interface::DBG_PRINT);

	// *********** WARNING ORDER IS IMPORTANT : LoadStateFromFile MUST BE CALLED FIRST! **************/

	// Loading the parameters
	ptrDH->loadStateFromFile(parametersFile);

	// Loading tree:
	std::vector<std::string> taxaNames = nexusParser->getTaxaNames();
	std::map<int, std::string> taxaIdToNameMap(nexusParser->getTaxaIdToNameMap());
	PS::TreeSharedPtr ptrTree( new PS::Tree(nexusParser) );
	std::string newickStr(ptrTree->getNewickString(taxaIdToNameMap));
	ptrDH->setTree(newickStr);

	// Trick to get probes automatically spaced
	if(ptrDH->synchMonitoring.size() == 1 && ptrDH->synchMonitoring.front() < 0.) {
		double step = fabs(ptrDH->synchMonitoring.front());
		ptrDH->synchMonitoring.clear();
		double ageOfOldest = ptrTree->getOldestNode()->getAge();

		for(double time=1.e-5; time < ageOfOldest; time += step) {
			ptrDH->synchMonitoring.push_back(time);
		}
		ptrDH->synchMonitoring.push_back(ageOfOldest-1.e-5);
	}

	// Loading data:
	if(nexusParser->hasData()) {
		Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

		// Loading the data
		std::map<std::string, std::vector<double> > probabilityMap;
		for(size_t iT=0; iT<taxaNames.size(); ++iT) {
			std::vector<double> prob(nexusParser->getNState(), 0.);
			// Load
			Eigen::VectorXd pVec = ptrData->getProbForTaxaLabel(taxaNames[iT]);
			// To std vector
			for(size_t iP=0; iP<prob.size(); ++iP) prob[iP] = pVec(iP);
			// in map
			probabilityMap[taxaNames[iT]] = prob;
		}

		ptrDH->setData(nexusParser->getTaxaNames(), probabilityMap);
	}

	Utils::Parallel::Manager::getInstance()->setMaxNThread(ptrDH->nThreads);
	Utils::Parallel::Manager::getInstance()->setNThread(ptrDH->nThreads);

	TensorPhylo::Interface::vecHistories_t vecHistories = ptrDH->drawMultipleHistories(nHistories);

	return vecHistories;

}

double computeLikForInterfaceTests(int idApproximator,
								   Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
								   const std::string &parametersFile) {

	//Utils::Parallel::Manager::getInstance()->setNThread(1);
	//Utils::Parallel::Manager::getInstance()->setMaxNThread(1);

	boost::shared_ptr<TensorPhylo::Interface::DistributionHandlerImpl> ptrDH(TensorPhylo::Interface::DistributionHandlerImpl::create());
	//ptrDH->setDebugMode(TensorPhylo::Interface::DBG_PRINT);

	// *********** WARNING ORDER IS IMPORTANT : LoadStateFromFile MUST BE CALLED FIRST! **************/

	// Loading the parameters
	ptrDH->loadStateFromFile(parametersFile);

	// Loading tree:
	std::vector<std::string> taxaNames = nexusParser->getTaxaNames();
	std::map<int, std::string> taxaIdToNameMap(nexusParser->getTaxaIdToNameMap());
	PS::TreeSharedPtr ptrTree( new PS::Tree(nexusParser) );
	std::string newickStr(ptrTree->getNewickString(taxaIdToNameMap));
	ptrDH->setTree(newickStr);

	assert(taxaNames.size() == ptrTree->getTerminalNodes().size() && "There are more taxa in the data file then the number of tips in the tree.");


	// Trick to get probes automatically spaced
	if(ptrDH->synchMonitoring.size() == 1 && ptrDH->synchMonitoring.front() < 0.) {
		double step = fabs(ptrDH->synchMonitoring.front());
		ptrDH->synchMonitoring.clear();
		double ageOfOldest = ptrTree->getOldestNode()->getAge();

		for(double time=1.e-5; time < ageOfOldest; time += step) {
			ptrDH->synchMonitoring.push_back(time);
		}
		ptrDH->synchMonitoring.push_back(ageOfOldest-1.e-5);
	}


	// Loading data:
	if(nexusParser->hasData()) {
		Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

		// Loading the data
		std::map<std::string, std::vector<double> > probabilityMap;
		for(size_t iT=0; iT<taxaNames.size(); ++iT) {
			std::vector<double> prob(nexusParser->getNState(), 0.);
			// Load
			Eigen::VectorXd pVec = ptrData->getProbForTaxaLabel(taxaNames[iT]);
			// To std vector
			for(size_t iP=0; iP<prob.size(); ++iP) prob[iP] = pVec(iP);
			// in map
			probabilityMap[taxaNames[iT]] = prob;
		}

		ptrDH->setData(nexusParser->getTaxaNames(), probabilityMap);
	} else { // If no data we init manually with ? chars.
		// Loading the data
		std::map<std::string, std::vector<double> > probabilityMap;
		for(size_t iT=0; iT<taxaNames.size(); ++iT) {
			std::vector<double> prob(ptrDH->rootPrior.size(), 1.);
			probabilityMap[taxaNames[iT]] = prob;
		}
		ptrDH->setData(taxaNames, probabilityMap);
	}

	// Ugly be necessary.
	TensorPhylo::Interface::approximatorVersion_t approxVersion = static_cast<TensorPhylo::Interface::approximatorVersion_t>(idApproximator);
	// Making sure that we test parelell approximator in pararalell

#ifdef _OPENMP
	bool changed = false;
	if(approxVersion == TensorPhylo::Interface::PARALLEL_BRANCHWISE || approxVersion == TensorPhylo::Interface::PARALLEL_OPTIMIZED) {
		if(ptrDH->nThreads == 1) {
			ptrDH->setNumberOfThreads(2);
			changed = true;
		}
		ptrDH->setLikelihoodApproximator(approxVersion);

		if(changed) {
			ptrDH->setNumberOfThreads(1);
		}

	} else {
		ptrDH->setLikelihoodApproximator(approxVersion);
	}
#else
	ptrDH->setLikelihoodApproximator(approxVersion);
#endif

	return ptrDH->computeLogLikelihood();

}

double computeLikForInterfaceTests(Likelihood::Integrator::integrationScheme_t intScheme,
								   const std::string &aNewickStringRB,
								   Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
								   Parameters::ContainerSharedPtr parameters) {

	boost::shared_ptr<TensorPhylo::Interface::DistributionHandlerImpl> ptrDH(TensorPhylo::Interface::DistributionHandlerImpl::create());

	// Loading tree:
	ptrDH->setTree(aNewickStringRB);

	// Loading data:
	Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

	std::vector<std::string> taxaNames(nexusParser->getTaxaNames());
	std::map<std::string, std::vector<double> > probabilityMap;
	for(size_t iT=0; iT<taxaNames.size(); ++iT) {
		std::vector<double> prob(nexusParser->getNState(), 0.);
		// Load
		Eigen::VectorXd pVec = ptrData->getProbForTaxaLabel(taxaNames[iT]);
		// To std vector
		for(size_t iP=0; iP<prob.size(); ++iP) prob[iP] = pVec(iP);
		// in map
		probabilityMap[taxaNames[iT]] = prob;
	}

	ptrDH->setData(nexusParser->getTaxaNames(), probabilityMap);


	{ // Loading misc:
		ptrDH->setInitialDeltaT(0.05);
		ptrDH->setApplyTreeLikCorrection(parameters->applyTreeLikCorrection);
		ptrDH->setLikelihoodApproximator(static_cast<TensorPhylo::Interface::approximatorVersion_t>(parameters->intLikApproximator));
		ptrDH->setConditionalProbabilityType(static_cast<TensorPhylo::Interface::conditionalProbability_t>(parameters->condType));
		ptrDH->setNumberOfThreads(parameters->nThreads);

		std::vector<double> rootPrior(parameters->rootPrior.size());
		for(size_t i=0; i<rootPrior.size(); ++i) rootPrior[i] = parameters->rootPrior(i);
		ptrDH->setRootPrior(rootPrior);
	}


	{ // Load parameters (vectors)
		std::vector<double> times;
		std::vector<double> lambda(parameters->lambda.data(), parameters->lambda.data()+parameters->lambda.size());
		TensorPhylo::Interface::stdMatrixXd matrixLambda(1, lambda);
		ptrDH->setLambda(times, matrixLambda);

		std::vector<double> mu(parameters->mu.data(), parameters->mu.data()+parameters->mu.size());
		TensorPhylo::Interface::stdMatrixXd matrixMu(1, mu);
		ptrDH->setMu(times, matrixMu);

		std::vector<double> phi(parameters->phi.data(), parameters->phi.data()+parameters->phi.size());
		TensorPhylo::Interface::stdMatrixXd matrixPhi(1, phi);
		ptrDH->setPhi(times, matrixPhi);

		std::vector<double> delta(parameters->delta.data(), parameters->delta.data()+parameters->delta.size());
		TensorPhylo::Interface::stdMatrixXd matrixDelta(1, delta);
		ptrDH->setDelta(times, matrixDelta);
	}

	{ // Load parameters (matrix, tensor)
		std::vector<double> times;
		TensorPhylo::Interface::stdMatrixXd eta;
		for(size_t i=0; i<(size_t)parameters->eta.rows(); ++i) {
			TensorPhylo::Interface::stdVectorXd rowEta;
			for(size_t j=0; j<(size_t)parameters->eta.cols(); ++j) {
				rowEta.push_back(parameters->eta(i,j));
			}
			eta.push_back(rowEta);
		}

		std::vector<TensorPhylo::Interface::stdMatrixXd> vecEta(1, eta);
		ptrDH->setEta(times, vecEta);

		TensorPhylo::Interface::eventMap_t eMap;
		for(size_t i=0; i<(size_t)parameters->omega.size(); ++i) {
			for(size_t j=0; j<(size_t)parameters->omega[i].rows(); ++j) {
				for(size_t k=0; k<(size_t)parameters->omega[i].cols(); ++k) {
					if(parameters->omega[i](j,k) != 0.) {
						std::vector<unsigned> key = boost::assign::list_of(i)(j)(k);
						eMap[key] = parameters->omega[i](j,k);
					}
				}
			}
		}

		std::vector<TensorPhylo::Interface::eventMap_t> vecOmega(1, eMap);
		ptrDH->setOmega(nexusParser->getNState(), times, vecOmega);
	}

	{ // Sync events: mass spec probs
		TensorPhylo::Interface::stdMatrixXd massSpecProbs;
		for(size_t i=0; i<parameters->massSpeciationProb.size(); ++i) {
			std::vector<double> vecProb(parameters->massSpeciationProb[i].data(), parameters->massSpeciationProb[i].data()+parameters->massSpeciationProb[i].size());
			massSpecProbs.push_back(vecProb);
		}
		ptrDH->setMassSpeciationEvents(parameters->massSpeciationTimes, massSpecProbs);
	}

	{ // Sync events: mass ext probs
		TensorPhylo::Interface::stdMatrixXd massExtProb;
		for(size_t i=0; i<parameters->massExtinctionProb.size(); ++i) {
			std::vector<double> vecProb(parameters->massExtinctionProb[i].data(), parameters->massExtinctionProb[i].data()+parameters->massExtinctionProb[i].size());
			massExtProb.push_back(vecProb);
		}
		ptrDH->setMassExtinctionEvents(parameters->massExtinctionTimes, massExtProb);

		std::vector< TensorPhylo::Interface::stdMatrixXd > massExtinctionStateChangeProb(parameters->massExtinctionStateChangeProb.size());
		for(size_t i=0; i<parameters->massExtinctionStateChangeProb.size(); ++i) {
			massExtinctionStateChangeProb[i].resize(parameters->massExtinctionStateChangeProb[i].rows());
			for(size_t j=0; j<(size_t)parameters->massExtinctionStateChangeProb[i].rows(); ++j) {
				massExtinctionStateChangeProb[i][j].resize(parameters->massExtinctionStateChangeProb[i].cols());
				for(size_t k=0; k<(size_t)parameters->massExtinctionStateChangeProb[i].cols(); ++k) {
					massExtinctionStateChangeProb[i][j][k] = parameters->massExtinctionStateChangeProb[i](j,k);
				}
			}
		}
		ptrDH->setMassExtinctionStateChangeProb(massExtinctionStateChangeProb);
	}

	{ // Sync events: mass sampling probs
		TensorPhylo::Interface::stdMatrixXd massSamplingProb;
		for(size_t i=0; i<parameters->massSamplingProb.size(); ++i) {
			std::vector<double> vecProb(parameters->massSamplingProb[i].data(), parameters->massSamplingProb[i].data()+parameters->massSamplingProb[i].size());
			massSamplingProb.push_back(vecProb);
		}
		ptrDH->setMassSamplingEvents(parameters->massSamplingTimes, massSamplingProb);
	}

	{ // Sync events: mass destr sampling probs
		TensorPhylo::Interface::stdMatrixXd massDestrSamplingProb;
		for(size_t i=0; i<parameters->massSpeciationProb.size(); ++i) {
			std::vector<double> vecProb(parameters->massDestrSamplingProb[i].data(), parameters->massDestrSamplingProb[i].data()+parameters->massDestrSamplingProb[i].size());
			massDestrSamplingProb.push_back(vecProb);
		}
		ptrDH->setMassDestrSamplingEvents(parameters->massDestrSamplingTimes, massDestrSamplingProb);
	}

	return ptrDH->computeLogLikelihood();
}

double computeLogLik(size_t nLikApproximation,
					 Likelihood::Integrator::integrationScheme_t intScheme,
					 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					 Parameters::ContainerSharedPtr parameters) {

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//Utils::Parallel::Manager::getInstance()->setNThread(1);
	//Utils::Parallel::Manager::getInstance()->setMaxNThread(1);

	PS::TreeSharedPtr ptrTree( new PS::Tree(nexusParser) );
	//ptrTree->scaleTree(5.);

	SynchronousEvents::ContainerSharedPtr syncEvents( new SynchronousEvents::Container(parameters) );

	// Create the tensors
	Tensor::ContainerSharedPtr tensors = Tensor::Factory::createContainerWithTimeHomogenousVectors(parameters);
	// We use the factory to create a tensor container with only time homogeneous vectors (based on the parameters)
	if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB)) {
		std::cout << "--------------------------" << std::endl << "Tensor properties : " << std::endl << tensors->reportTensorInformations()  << "--------------------------" << std::endl;
	}

	// Trick to get probes automatically spaced
	if(parameters->synchMonitoring.size() == 1 && parameters->synchMonitoring.front() < 0.) {
		double step = fabs(parameters->synchMonitoring.front());
		parameters->synchMonitoring.clear();
		double ageOfOldest = ptrTree->getOldestNode()->getAge();

		for(double time=1.e-5; time < ageOfOldest; time += step) {
			parameters->synchMonitoring.push_back(time);
		}
		parameters->synchMonitoring.push_back(ageOfOldest-1.e-5);
	}

	assert((size_t)parameters->lambda.size() == nexusParser->getNState() && "The number of state indicated in the parameters file does not aggree with the data.");
	Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

	// Init scheduler
	using Likelihood::Scheduler::BaseScheduler;
	boost::shared_ptr<BaseScheduler> scheduler( new BaseScheduler(ptrTree) );
	// Set synchronous events in the scheduler
	scheduler->setSynchronousEvents(syncEvents);
	// Set probes
	scheduler->setMonitoringProbes(parameters->synchMonitoring);
	scheduler->defineAndSetRescalingEvents();

	/*for(size_t i=0; i<scheduler->getEvents().size(); ++i) {
		std::cout << scheduler->getEvents()[i]->toString() << std::endl;
	}*/

	// Vague prior
	Eigen::VectorXd prior = parameters->rootPrior;

	// Create the likelihood approximator
	Likelihood::Conditions::conditionalProbability_t condType = Likelihood::Conditions::intToConditionalProbabilityType(parameters->condType);

	assert(parameters->intLikApproximator >= 0 && (size_t)parameters->intLikApproximator < Likelihood::Approximator::APPROXIMATOR_NAMES.size());
	Likelihood::Approximator::approximatorVersion_t approximatorVersion = static_cast<Likelihood::Approximator::approximatorVersion_t>(parameters->intLikApproximator);
	Likelihood::Approximator::ApproximatorSharedPtr approximator;
	if(approximatorVersion == Likelihood::Approximator::SEQUENTIAL_OPTIMIZED) {
		approximator = (Likelihood::Approximator::Factory::createSequentialTemplateCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
	} else if(approximatorVersion == Likelihood::Approximator::SEQUENTIAL_BRANCHWISE) {
		approximator = (Likelihood::Approximator::Factory::createSequentialBranchwiseCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
#if defined(_OPENMP)
	} else if(approximatorVersion == Likelihood::Approximator::PARALLEL_OPTIMIZED) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(parameters->nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(parameters->nThreads);
		//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
		approximator = (Likelihood::Approximator::Factory::createParallelOpenMPCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
	} else if(approximatorVersion == Likelihood::Approximator::PARALLEL_BRANCHWISE) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(parameters->nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(parameters->nThreads);
		//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
		approximator = (Likelihood::Approximator::Factory::createParallelBranchwiseCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
#endif
	} else if(approximatorVersion == Likelihood::Approximator::AUTO_TUNING) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(parameters->nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(parameters->nThreads);
		approximator = (Likelihood::Approximator::Factory::createAutoTuningCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
	} else {
		assert(approximator != NULL && "The approximator requested is not available (openmp approximators are only available if you compiled with openmp (-fopenmp)).");
	}

	if(parameters->applyTreeLikCorrection) {
		approximator->enableTreeLikelihoodCorrection();
	} else {
		approximator->disableTreeLikelihoodCorrection();
	}

	approximator->setDefaultDeltaT(parameters->deltaT);
	// If not done will crash
	approximator->setPriorStateProbability(prior);

	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	if(Utils::Output::outputManager().check(Utils::Output::LOW_VERB)) {
		std::cout << "Computing " << nLikApproximation << " likelihoods." << std::endl;
	}

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")

	double logLik = 0.;
	if(nLikApproximation == 1) {
		logLik = approximator->approximateLogLikelihood();
	} else {
		for(size_t iL=0; iL<nLikApproximation; ++iL) {
			_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "LIK_TIME")
			logLik = approximator->approximateLogLikelihood();
			_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "LIK_TIME")
			if((nLikApproximation > 40) && (iL % (size_t)(nLikApproximation/20.0)) == 0) {
				std::cout << "*";
				std::cout.flush();
			}
		}
		std::cout << std::endl;
	}


	_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")


	if(Utils::Output::outputManager().check(Utils::Output::LOW_VERB)) {
		std::cout << "LogLik = " << logLik << std::scientific << std::setprecision(8) << " ( scientific: " << logLik << " ) " <<  std::endl;
		std::cout << std::fixed << std::setprecision(6) <<  std::endl;
		std::cout << "Total number of integration steps : " << approximator->getTotalNumberOfIntegrationSteps() << std::endl;
	}

	if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB)) {
		approximator->orderProbes();
		const std::vector<Likelihood::Monitor::ProbeState>& probes = approximator->getObservedProbesState();
		std::cout << "------- PROBES --------" << std::endl;
		for(size_t i=0; i<probes.size(); ++i) {
			std::cout << probes[i].toString() << std::endl;
		}
	}

	if(Utils::Output::outputManager().check(Utils::Output::DUMP_VERB)) {
		const std::vector<Likelihood::Monitor::ProbeState>& probes = approximator->getObservedProbesState();
		std::cout << "------- DUMPING PROBES IN 'dump_probes.txt' --------" << std::endl;
		std::ofstream oFile(Utils::Output::outputManager().getDumpFileName());
		for(size_t i=0; i<probes.size(); ++i) {
			oFile << probes[i].toDumpFormat();
		}
	}

	_REPORT_TO_FILE(Utils::Profiling::UniqueProfiler::getInstance())

	return logLik;
}

TensorPhylo::Interface::vecHistories_t computeStochasticMappingLegacy(size_t nReplicas,
					 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					 Parameters::ContainerSharedPtr parameters,
					 PS::TreeSharedPtr ptrTree) {
	return computeStochasticMappingLegacy(nReplicas, Likelihood::Approximator::StochasticMapping::REJECTION_SAMPLING_ALGO, nexusParser, parameters, ptrTree);
}

TensorPhylo::Interface::vecHistories_t computeStochasticMappingLegacy(size_t nReplicas,
		 	 	 	 Likelihood::Approximator::StochasticMapping::stochasticMappingAlgo_t aAlgo,
					 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					 Parameters::ContainerSharedPtr parameters,
					 PS::TreeSharedPtr ptrTree) {

	_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//Utils::Parallel::Manager::getInstance()->setNThread(1);
	//Utils::Parallel::Manager::getInstance()->setMaxNThread(1);

	//ptrTree->scaleTree(5.);

	SynchronousEvents::ContainerSharedPtr syncEvents( new SynchronousEvents::Container(parameters) );

	// Create the tensors
	Tensor::ContainerSharedPtr tensors = Tensor::Factory::createContainerWithTimeHomogenousVectors(parameters);
	// We use the factory to create a tensor container with only time homogeneous vectors (based on the parameters)
	if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB)) {
		std::cout << "--------------------------" << std::endl << "Tensor properties : " << std::endl << tensors->reportTensorInformations()  << "--------------------------" << std::endl;
	}

	// Trick to get probes automatically spaced
	if(parameters->synchMonitoring.size() == 1 && parameters->synchMonitoring.front() < 0.) {
		double step = fabs(parameters->synchMonitoring.front());
		parameters->synchMonitoring.clear();
		double ageOfOldest = ptrTree->getOldestNode()->getAge();

		for(double time=1.e-5; time < ageOfOldest; time += step) {
			parameters->synchMonitoring.push_back(time);
		}
		parameters->synchMonitoring.push_back(ageOfOldest-1.e-5);
	}

	assert((size_t)parameters->lambda.size() == nexusParser->getNState() && "The number of state indicated in the parameters file does not aggree with the data.");
	Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

	// Init scheduler
	using Likelihood::Scheduler::BaseScheduler;
	boost::shared_ptr<BaseScheduler> scheduler( new BaseScheduler(ptrTree) );
	// Set synchronous events in the scheduler
	scheduler->setSynchronousEvents(syncEvents);
	// Set probes
	scheduler->setMonitoringProbes(parameters->synchMonitoring);
	scheduler->defineAndSetRescalingEvents();

	/*for(size_t i=0; i<scheduler->getEvents().size(); ++i) {
		std::cout << scheduler->getEvents()[i]->toString() << std::endl;
	}*/

	// Vague prior
	Eigen::VectorXd prior = parameters->rootPrior;

	// Create the likelihood approximator
	Likelihood::Conditions::conditionalProbability_t condType = Likelihood::Conditions::intToConditionalProbabilityType(parameters->condType);

	Likelihood::Approximator::StochasticMapping::StochasticMappingApproxSharedPtr ptrStochasticMapper =
			Likelihood::Approximator::Factory::createStochasticMappingApprox(Likelihood::Integrator::RUNGE_KUTTA_DOPRI5, condType, ptrData, scheduler, syncEvents, tensors);

	ptrStochasticMapper->setAlgorithm(aAlgo);

	if(parameters->applyTreeLikCorrection) {
		ptrStochasticMapper->enableTreeLikelihoodCorrection();
	} else {
		ptrStochasticMapper->disableTreeLikelihoodCorrection();
	}

	ptrStochasticMapper->setDefaultDeltaT(parameters->deltaT);
	// If not done will crash
	ptrStochasticMapper->setPriorStateProbability(prior);

	ptrStochasticMapper->approximateLogLikelihood();

	TensorPhylo::Interface::vecHistories_t vecHistories;
	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	for(size_t iR=0; iR<nReplicas; ++iR) {
		vecHistories.push_back(ptrStochasticMapper->drawHistory());
	}


	if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB)) {
		for(TensorPhylo::Interface::mapHistories_t::iterator itH = vecHistories.front().begin(); itH != vecHistories.front().end(); ++itH) {

			std::cout << "History for edge " << itH->first << " : " << std::endl;
			for(size_t iT=0; iT<itH->second.size(); ++iT) {
				std::cout << itH->second[iT].first << " -> " << itH->second[iT].second << std::endl;
			}
			std::cout << "--------------------------------------" << std::endl;
		}
	}

	return vecHistories;

}


void runBenchmarkLegacy(size_t nLikApproximation,
					 Likelihood::Integrator::integrationScheme_t intScheme,
					 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					 Parameters::ContainerSharedPtr parameters,
					 const std::string &logFile) {

	std::stringstream ssReport;

	Utils::Profiling::CustomProfiling cp;

	PS::TreeSharedPtr ptrTree( new PS::Tree(nexusParser) );
	//ptrTree->scaleTree(5.);

	SynchronousEvents::ContainerSharedPtr syncEvents( new SynchronousEvents::Container(parameters) );

	// Create the tensors
	Tensor::ContainerSharedPtr tensors = Tensor::Factory::createContainerWithTimeHomogenousVectors(parameters);
	// We use the factory to create a tensor container with only time homogeneous vectors (based on the parameters)

	assert((size_t)parameters->lambda.size() == nexusParser->getNState() && "The number of state indicated in the parameters file does not aggree with the data.");
	Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

	// Init scheduler
	using Likelihood::Scheduler::BaseScheduler;
	boost::shared_ptr<BaseScheduler> scheduler( new BaseScheduler(ptrTree) );
	// Set synchronous events in the scheduler
	scheduler->setSynchronousEvents(syncEvents);
	// Set probes
	scheduler->setMonitoringProbes(parameters->synchMonitoring);
	scheduler->defineAndSetRescalingEvents();

	// Vague prior
	Eigen::VectorXd prior = parameters->rootPrior;

	// Create the likelihood approximator
	Likelihood::Conditions::conditionalProbability_t condType = Likelihood::Conditions::intToConditionalProbabilityType(parameters->condType);

	std::cout << "Benchmarking [" << Likelihood::Approximator::APPROXIMATOR_NAMES[parameters->intLikApproximator] << "]" << std::endl;

	cp.startTime(0);

	assert(parameters->intLikApproximator >= 0 && (size_t)parameters->intLikApproximator < Likelihood::Approximator::APPROXIMATOR_NAMES.size());
	Likelihood::Approximator::approximatorVersion_t approximatorVersion = static_cast<Likelihood::Approximator::approximatorVersion_t>(parameters->intLikApproximator);
	Likelihood::Approximator::ApproximatorSharedPtr approximator;

	if(approximatorVersion == Likelihood::Approximator::SEQUENTIAL_OPTIMIZED) {
		approximator = (Likelihood::Approximator::Factory::createSequentialTemplateCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
	} else if(approximatorVersion == Likelihood::Approximator::SEQUENTIAL_BRANCHWISE) {
		approximator = (Likelihood::Approximator::Factory::createSequentialBranchwiseCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
#if defined(_OPENMP)
	} else if(approximatorVersion == Likelihood::Approximator::PARALLEL_OPTIMIZED) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(parameters->nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(parameters->nThreads);
		approximator = (Likelihood::Approximator::Factory::createParallelOpenMPCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
	} else if(approximatorVersion == Likelihood::Approximator::PARALLEL_BRANCHWISE) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(parameters->nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(parameters->nThreads);
		//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
		approximator = (Likelihood::Approximator::Factory::createParallelBranchwiseCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
#endif
	} else if(approximatorVersion == Likelihood::Approximator::AUTO_TUNING) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(parameters->nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(parameters->nThreads);
		approximator = (Likelihood::Approximator::Factory::createAutoTuningCPU(intScheme, condType, ptrData, scheduler, syncEvents, tensors));
	} else {
		assert(approximator != NULL && "The approximator requested is not available (openmp approximators are only available if you compiled with openmp (-fopenmp)).");
	}

	if(parameters->applyTreeLikCorrection) {
		approximator->enableTreeLikelihoodCorrection();
	} else {
		approximator->disableTreeLikelihoodCorrection();
	}

	approximator->setDefaultDeltaT(parameters->deltaT);
	// If not done will crash
	approximator->setPriorStateProbability(prior);

	cp.endTime(0);

	double logLik = 0.;
	cp.startTime(1);
	for(size_t iL=0; iL<nLikApproximation; ++iL) {
		logLik = approximator->approximateLogLikelihood();
	}
	cp.endTime(1);
	ssReport << std::scientific << Likelihood::Approximator::APPROXIMATOR_NAMES[parameters->intLikApproximator] << "\t" << 1 << "\t" << cp.duration(0) << "\t" << cp.duration(1);


	ssReport << "\t" << std::scientific << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << logLik << std::endl;

	// Write result.
	std::ofstream oFile(logFile.c_str(), std::ios::out);
	assert(oFile.good());
	oFile << ssReport.str() << std::endl;
	oFile.close();

}

void runBenchmarkNew(size_t nLikApproximation,
					 Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
					 const std::string &parametersFile,
					 const std::string &logFile) {

	std::stringstream ssReport;

	Utils::Profiling::CustomProfiling cp;


	cp.startTime(0);

	boost::shared_ptr<TensorPhylo::Interface::DistributionHandlerImpl> ptrDH(TensorPhylo::Interface::DistributionHandlerImpl::create());
	//ptrDH->setDebugMode(TensorPhylo::Interface::DBG_PRINT);

	// *********** WARNING ORDER IS IMPORTANT : LoadStateFromFile MUST BE CALLED FIRST! **************/

	// Loading the parameters
	ptrDH->loadStateFromFile(parametersFile);

	// Loading tree:
	std::vector<std::string> taxaNames = nexusParser->getTaxaNames();
	std::map<int, std::string> taxaIdToNameMap(nexusParser->getTaxaIdToNameMap());
	PS::TreeSharedPtr ptrTree( new PS::Tree(nexusParser) );
	std::string newickStr(ptrTree->getNewickString(taxaIdToNameMap));
	ptrDH->setTree(newickStr);


	// Loading data:
	if(nexusParser->hasData()) {
		Phylogeny::Data::ContainerSharedPtr ptrData = nexusParser->createDataContainer();

		// Loading the data
		std::map<std::string, std::vector<double> > probabilityMap;
		for(size_t iT=0; iT<taxaNames.size(); ++iT) {
			std::vector<double> prob(nexusParser->getNState(), 0.);
			// Load
			Eigen::VectorXd pVec = ptrData->getProbForTaxaLabel(taxaNames[iT]);
			// To std vector
			for(size_t iP=0; iP<prob.size(); ++iP) prob[iP] = pVec(iP);
			// in map
			probabilityMap[taxaNames[iT]] = prob;
		}

		ptrDH->setData(nexusParser->getTaxaNames(), probabilityMap);
	}

	cp.endTime(0);

	Utils::Parallel::Manager::getInstance()->setMaxNThread(ptrDH->nThreads);
	Utils::Parallel::Manager::getInstance()->setNThread(ptrDH->nThreads);

	double logLik = 0.;
	std::cout << "Benchmarking." << std::endl;

	cp.startTime(1);
	for(size_t iL=0; iL<nLikApproximation; ++iL) {
		logLik = ptrDH->computeLogLikelihood();
	}
	cp.endTime(1);
	ssReport << std::scientific << "TensorPhyloViaInterface" << "\t" << cp.duration(0) << "\t" << cp.duration(1);

	ssReport << "\t" << std::scientific << logLik << "\t" << ptrDH->approximator->getTotalNumberOfIntegrationSteps() << std::endl;

	// Write result.
	std::ofstream oFile(logFile.c_str(), std::ios::out);
	assert(oFile.good());
	oFile << ssReport.str() << std::endl;
	oFile.close();

}


void benchmarkConvolution(const size_t M, const std::vector<size_t> &vecK) {

	Utils::Output::outputManager().setVerbosityThreshold(Utils::Output::SILENT_VERB);

	for(std::vector<size_t>::const_iterator it = vecK.begin(); it != vecK.end(); ++it) {

		size_t K = *it;

		std::stringstream ssK;
		ssK << "Testing 'band' matrix - vector multiplication for size K = " << K;
		std::cout << ssK.str() << std::endl;

		// Init the dense matrix
		double alpha = 42.5;
		Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(K,K);

		for(size_t i=0; i<K; ++i) {
			double sum = 0.;

			if(i > 0) {
				dense(i,i-1) = alpha;
				sum += alpha;
			}

			if(i < K-1) {
				dense(i,i+1) = alpha;
				sum += alpha;
			}

			dense(i,i) = -sum;
		}

		// Init the sparse matrix
		Tensor::eigenSparseMatrix_t sparseMat = dense.sparseView();
		sparseMat.makeCompressed();

		// Init vectors
		Eigen::VectorXd vecOnes = Eigen::VectorXd::Ones(K);
		Eigen::VectorXd vecRand = Eigen::VectorXd::Random(K);

		Utils::Profiling::CustomProfiling cp;


		cp.startTime();
		Eigen::VectorXd refOnes, refRand;
		for(size_t i=0; i<M; ++i) {
			refOnes = dense * vecOnes;
			refRand = dense * vecRand;
		}
		cp.endTime();
		double timeDense = cp.duration();
		std::cout << "Dense Eigen : t = " << timeDense << std::endl;



		{
			Eigen::VectorXd resOnes, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resOnes = sparseMat * vecOnes;
				resRand = sparseMat * vecRand;
			}
			cp.endTime();

			assert((refOnes-resOnes).norm() < 1.e-7);
			assert((refRand-resRand).norm() < 1.e-7);

		}
		std::cout << "Sparse Eigen : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration() << std::endl;

		{
			Eigen::VectorXd resOnes, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resOnes = Likelihood::Kernels::CPU::computeQuasseEtaVector<1>(alpha, vecOnes);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaVector<1>(alpha, vecRand);
			}
			cp.endTime();

			assert((refOnes-resOnes).norm() < 1.e-7);
			assert((refRand-resRand).norm() < 1.e-7);
		}
		std::cout << "Convolution stride = 1 : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::VectorXd resOnes, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resOnes = Likelihood::Kernels::CPU::computeQuasseEtaVector<2>(alpha, vecOnes);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaVector<2>(alpha, vecRand);
			}
			cp.endTime();

			assert((refOnes-resOnes).norm() < 1.e-7);
			assert((refRand-resRand).norm() < 1.e-7);
		}
		std::cout << "Convolution stride = 2 : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::VectorXd resOnes, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resOnes = Likelihood::Kernels::CPU::computeQuasseEtaVector<4>(alpha, vecOnes);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaVector<4>(alpha, vecRand);
			}
			cp.endTime();

			assert((refOnes-resOnes).norm() < 1.e-7);
			assert((refRand-resRand).norm() < 1.e-7);
		}
		std::cout << "Convolution stride = 4 : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::VectorXd resOnes, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resOnes = Likelihood::Kernels::CPU::computeQuasseEtaVector<6>(alpha, vecOnes);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaVector<6>(alpha, vecRand);
			}
			cp.endTime();

			assert((refOnes-resOnes).norm() < 1.e-7);
			assert((refRand-resRand).norm() < 1.e-7);
		}
		std::cout << "Convolution stride = 6 : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::VectorXd resOnes, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resOnes = Likelihood::Kernels::CPU::computeQuasseEtaVector<8>(alpha, vecOnes);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaVector<8>(alpha, vecRand);
			}
			cp.endTime();

			assert((refOnes-resOnes).norm() < 1.e-7);
			assert((refRand-resRand).norm() < 1.e-7);
		}
		std::cout << "Convolution stride = 8 : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;
		std::cout << "----------------------------------------------------------" << std::endl;
	}
}

void benchmarkConvolution2(const size_t M, const std::vector<size_t> &vecK) {

	Utils::Output::outputManager().setVerbosityThreshold(Utils::Output::SILENT_VERB);

	for(std::vector<size_t>::const_iterator it = vecK.begin(); it != vecK.end(); ++it) {

		size_t K = *it;

		std::stringstream ssK;
		ssK << "Testing 'band' matrix - matrix multiplication for size K = " << K;
		std::cout << ssK.str() << std::endl;

		// Init the dense matrix
		double alpha = 42.5;
		Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(K,K);

		for(size_t i=0; i<K; ++i) {
			double sum = 0.;

			if(i > 0) {
				dense(i,i-1) = alpha;
				sum += alpha;
			}

			if(i < K-1) {
				dense(i,i+1) = alpha;
				sum += alpha;
			}

			dense(i,i) = -sum;
		}

		// Init the sparse matrix
		Tensor::eigenSparseMatrix_t sparseMat = dense.sparseView();
		sparseMat.makeCompressed();

		// Init vectors
		Eigen::MatrixXd matI = Eigen::MatrixXd::Identity(K, K);
		Eigen::MatrixXd matRand = Eigen::MatrixXd::Random(K, K);


		Utils::Profiling::CustomProfiling cp;


		cp.startTime();
		Eigen::MatrixXd refI, refRand;
		for(size_t i=0; i<M; ++i) {
			refI = dense * matI;
			refRand = dense * matRand;
		}
		cp.endTime();
		double timeDense = cp.duration();
		std::cout << "Dense Eigen : t = " << timeDense << std::endl;

		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = sparseMat * matI;
				resRand = sparseMat * matRand;
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);

		}
		std::cout << "Sparse Eigen : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration() << std::endl;


		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<4, 1>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<4, 1>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 4 x 1  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<4, 4>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<4, 4>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 4 x 4  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<4, 8>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<4, 8>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 4 x 8  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;


		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<8, 4>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<8, 4>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 8 x 4  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<8, 8>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<8, 8>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 8 x 8  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;


		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<16, 4>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<16, 4>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 16 x 4  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<16, 8>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<16, 8>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 16 x 8  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;

		{
			Eigen::MatrixXd resI, resRand;
			cp.startTime();
			for(size_t i=0; i<M; ++i) {
				resI = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<16, 16>(alpha, matI);
				resRand = Likelihood::Kernels::CPU::computeQuasseEtaMatrix<16, 16>(alpha, matRand);
			}
			cp.endTime();

			assert((resI-refI).colwise().norm().sum() < 1.e-7);
			assert((resRand-refRand).colwise().norm().sum() < 1.e-7);
		}
		std::cout << "Convolution stride = 16 x 16  : t = " << cp.duration() << " :: speedup = " << timeDense/cp.duration()  << std::endl;


		std::cout << "----------------------------------------------------------" << std::endl;
	}
}

bool checkHistoriesEndStatesAgainstDataProbability(
		PS::TreeSharedPtr ptrTree,
		Phylogeny::NexusReader::NexusParserSharedPtr nexusParser,
		TensorPhylo::Interface::vecHistories_t vecHistories,
		double tolerance) {

	bool areCorrect = true;

	const std::vector<PS::Node*>& nodes = ptrTree->getNodes();
	for(size_t iN=0; iN<nodes.size(); ++iN) {								// For each Node
		int iTaxa = nodes[iN]->getTaxaId();
		if(iTaxa < 0) continue; 											// For each Node having data

		// Estimate observaed state frequency
		Eigen::VectorXd observedStateFreq =
				Eigen::VectorXd::Zero(nexusParser->getNState());
		for(size_t iV=0; iV<vecHistories.size(); ++iV) {					// For each sample history
			TensorPhylo::Interface::mapHistories_t::iterator itH;
			itH = vecHistories[iV].find(nodes[iN]->getId());
			assert(itH != vecHistories[iV].end());
			size_t state = itH->second.back().second;						// Last state for the parent branch

			observedStateFreq(state) += 1;
		}
		observedStateFreq /= observedStateFreq.sum();

		// Compute error with data probability
		double error = 0.;
		const std::vector<double> &dataStateFreq = nexusParser->getProbForTaxaId(iTaxa);
		double dataStateFreqSum = 0.;
		for(size_t iS=0; iS<dataStateFreq.size(); ++iS) {
			dataStateFreqSum += dataStateFreq[iS];
		}
		for(size_t iS=0; iS<dataStateFreq.size(); ++iS) {
			error += fabs(dataStateFreq[iS]/dataStateFreqSum - observedStateFreq(iS));
		}
		if(error > tolerance) {
			std::cout << "iTaxa = " << iTaxa << std::endl;
			std::cout << "observedStateFreq = " << observedStateFreq.transpose() << std::endl;
			std::cout << "dataStateFreq = ";
			for(size_t iS=0; iS<dataStateFreq.size(); ++iS) {
				std::cout << dataStateFreq[iS] << ", ";
			}
			std::cout << std::endl;
			return false;
		}
	}

	return areCorrect;
}

} /* namespace SequentialCPU */
} /* namespace Test */
