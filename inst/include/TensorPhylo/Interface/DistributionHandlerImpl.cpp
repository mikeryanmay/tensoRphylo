/*
 * DistributionHandlerImpl.cpp
 *
 *  Created on: Mar 9, 2020
 *      Author: xaviermeyer
 */

#include "DistributionHandlerImpl.h"

#include "Utils/ParametersSerializer.h"
#include "Utils/Parallel/Manager.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Parameters/IncParameterContainer.h"
#include "Data/Reader/IncPhyloReader.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Likelihood/ConditionTypes/ConditionType.h"
#include "Likelihood/Approximator/IncLikelihoodApproximator.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "../Likelihood/Approximator/StochasticMapping/SequentialStochasticMappingCPU.h"
#include "../Utils/RNG/Manager.h"
#include "Tensor/IncTensor.h"
#include "Likelihood/Kernels/CPU/Misc/QuasistationaryFrequency.h"

#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION
#include <boost/dll/alias.hpp>

BOOST_DLL_ALIAS(
	TensorPhylo::Interface::DistributionHandlerImpl::create,
    createTensorPhyloDistributionHandler);


namespace TensorPhylo {
namespace Interface {

DistributionHandlerImpl::DistributionHandlerImpl() :
		DistributionHandler(),
		debugMode(DBG_NONE), schedulerOperation(NONE),
		nThreads(1), initDeltaT(0.05),
		applyTreeLikCorrection(true), useQuasistationaryFrequency(false), compatibilityMode(false),
		approxVersion(approximatorVersion_t::AUTO_TUNING),
		condProbType(conditionalProbability_t::TIME),
		integrationScheme(integrationScheme_t::RUNGE_KUTTA_DOPRI5),
		ptrData(new Phylogeny::Data::Container()),
		ptrAsyncParams(new Parameters::AsyncParameterContainer()),
		ptrSyncParams(new Parameters::SyncParameterContainer()) {

	dirtyAsyncRateShifts = dirtySyncEventTimes = dirtyApproximator = false;
	dirtyLambda = dirtyMu = dirtyDelta = dirtyPhi = dirtyEta = dirtyOmega = false;

}

DistributionHandlerImpl::~DistributionHandlerImpl() {
}

void DistributionHandlerImpl::setTree(const std::string &aNewickTreeStr) {

	Phylogeny::NewickReader::NewickParserSharedPtr ptrNewickParser(new Phylogeny::NewickReader::NewickParser(aNewickTreeStr, Phylogeny::NewickReader::NewickParser::IS_TREE_STRING));
	ptrTree.reset( new Phylogeny::Structure::Tree(ptrNewickParser) );

	ptrData->setTaxaToIdMap(ptrNewickParser->getTaxaNames());

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] RECEIVED TREE: " << aNewickTreeStr << std::endl;
		ss << "[DEBUG TENSORPHYLO] TRANSLATED TO TREE: " << ptrTree->getNewickString() << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

	schedulerOperation = RESET;
}

void DistributionHandlerImpl::setData(const std::vector<std::string> &aTaxa, const std::map<std::string, std::vector<double> > &aProbabilityMap) {

	for(size_t iT = 0; iT < aTaxa.size(); ++iT ) {
		std::map<std::string, std::vector<double> >::const_iterator itFind = aProbabilityMap.find(aTaxa[iT]);
		assert(itFind != aProbabilityMap.end());
		ptrData->registerTaxaData(aTaxa[iT], itFind->second);
	}

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] RECEIVED DATA: " << std::endl;
		for(size_t iT = 0; iT < aTaxa.size(); ++iT ) {
			ss << "TAXA :" << aTaxa[iT] << std::endl;
			std::map<std::string, std::vector<double> >::const_iterator itFind = aProbabilityMap.find(aTaxa[iT]);
			assert(itFind != aProbabilityMap.end());
			ss << "RB DATA : ";
			for(size_t i=0; i<itFind->second.size(); ++i) {
				ss << itFind->second[i] <<" ";
			}
			ss << std::endl;
			ss << "TP DATA : " << ptrData->getProbForTaxaLabel(aTaxa[iT]).transpose() << std::endl;
		}
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());

	}

	schedulerOperation = RESET;
}

void DistributionHandlerImpl::forceSchedulerUpdate() {
	schedulerOperation = UPDATE;
}

void DistributionHandlerImpl::forceApproximatorDirty() {
	dirtyApproximator = true;
}


void DistributionHandlerImpl::setApplyTreeLikCorrection(bool aDoApply) {
	applyTreeLikCorrection = aDoApply;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] ";
		if(applyTreeLikCorrection) ss << "Applying tree lik correction!" << std::endl;
		else ss << "NOT applying tree lik correction!" << std::endl;
		ss << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setConditionalProbCompatibilityMode(bool setActive) {
	compatibilityMode = setActive;
	dirtyApproximator = true;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] ";
		if(compatibilityMode) ss << "Swapping to compatibility mode!" << std::endl;
		else ss << "Disabling compatibility mode!" << std::endl;
		ss << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}
}

void DistributionHandlerImpl::setQuasistationaryFrequencyMode(bool setActive) {
	useQuasistationaryFrequency = setActive;
	dirtyApproximator = true;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] ";
		if(useQuasistationaryFrequency) ss << "Using quasistationary frequency at the root!" << std::endl;
		else ss << "NOT using quasistationary frequency at the root!" << std::endl;
		ss << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setLikelihoodApproximator(approximatorVersion_t aApproxVersion) {
	approxVersion = aApproxVersion;
	dirtyApproximator = true;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] Approximator version = " << static_cast<int>(approxVersion) << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setConditionalProbabilityType(conditionalProbability_t aCondProb) {
	condProbType = aCondProb;
	dirtyApproximator = true;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] Conditional probability type = " << static_cast<int>(condProbType) << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setIntegrationScheme(integrationScheme_t aIntScheme) {
	integrationScheme = aIntScheme;
	dirtyApproximator = true;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] Integration scheme type = " << static_cast<int>(integrationScheme) << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setNumberOfThreads(size_t aNThreads) {
	nThreads = aNThreads;

	Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
	Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);
	dirtyApproximator = true;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] N Threads  = " << nThreads << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}
}

void DistributionHandlerImpl::setInitialDeltaT(double aInitDeltaT) {
	initDeltaT = aInitDeltaT;

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] initDeltaT  = " << aInitDeltaT << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setRootPrior(const stdVectorXd &aRootPrior) {

	if(rootPrior.size() == 0) {
		rootPrior.resize(aRootPrior.size());
	}

	for(size_t iP=0; iP<aRootPrior.size(); ++iP) {
		rootPrior(iP) = aRootPrior[iP];
	}

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] ROOT FREQ: " << rootPrior.transpose() << std::endl;
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setLambda(const stdVectorXd &aTimes, const stdMatrixXd &aLambda) {

	bool areTimesDifferent = copyVectorWithCheck(aTimes, ptrAsyncParams->timesLambda);
	bool areRatesDifferent = copyMatrixWithCheck(aLambda, ptrAsyncParams->vecLambda);

	if(areTimesDifferent) {
		dirtyAsyncRateShifts = true;
		updateSchedulerOperationTimeChange();
	}
	if(areRatesDifferent) {
		dirtyLambda = true;
		updateSchedulerOperationRateChange();
	}

	printParamVectorDebug("Lambda", ptrAsyncParams->timesLambda, ptrAsyncParams->vecLambda);

}

void DistributionHandlerImpl::setMu(const stdVectorXd &aTimes, const stdMatrixXd &aMu) {

	bool areTimesDifferent = copyVectorWithCheck(aTimes, ptrAsyncParams->timesMu);
	bool areRatesDifferent = copyMatrixWithCheck(aMu, ptrAsyncParams->vecMu);

	if(areTimesDifferent) {
		dirtyAsyncRateShifts = true;
		updateSchedulerOperationTimeChange();
	}
	if(areRatesDifferent) {
		dirtyMu = true;
		updateSchedulerOperationRateChange();
	}

	printParamVectorDebug("Mu", ptrAsyncParams->timesMu, ptrAsyncParams->vecMu);

}

void DistributionHandlerImpl::setPhi(const stdVectorXd &aTimes, const stdMatrixXd &aPhi) {

	bool areTimesDifferent = copyVectorWithCheck(aTimes, ptrAsyncParams->timesPhi);
	bool areRatesDifferent = copyMatrixWithCheck(aPhi, ptrAsyncParams->vecPhi);

	if(areTimesDifferent) {
		dirtyAsyncRateShifts = true;
		updateSchedulerOperationTimeChange();
	}
	if(areRatesDifferent) {
		dirtyPhi = true;
		updateSchedulerOperationRateChange();
	}

	printParamVectorDebug("Phi", ptrAsyncParams->timesPhi, ptrAsyncParams->vecPhi);

}

void DistributionHandlerImpl::setDelta(const stdVectorXd &aTimes, const stdMatrixXd &aDelta) {

	bool areTimesDifferent = copyVectorWithCheck(aTimes, ptrAsyncParams->timesDelta);
	bool areRatesDifferent = copyMatrixWithCheck(aDelta, ptrAsyncParams->vecDelta);

	if(areTimesDifferent) {
		dirtyAsyncRateShifts = true;
		updateSchedulerOperationTimeChange();
	}
	if(areRatesDifferent) {
		dirtyDelta = true;
		updateSchedulerOperationRateChange();
	}

	printParamVectorDebug("Delta", ptrAsyncParams->timesDelta, ptrAsyncParams->vecDelta);

}

void DistributionHandlerImpl::setEta(const stdVectorXd &aTimes, const std::vector< stdMatrixXd > &aEta) {

	bool areTimesDifferent = copyVectorWithCheck(aTimes, ptrAsyncParams->timesEta);
	bool areRatesDifferent = copyVectorOfRateMatrixWithCheck(aEta, ptrAsyncParams->vecEta);

	if(areTimesDifferent) {
		dirtyAsyncRateShifts = true;
		updateSchedulerOperationTimeChange();
	}
	if(areRatesDifferent) {
		dirtyEta = true;
		updateSchedulerOperationRateChange();
	}

	printParamMatrixDebug("Eta", ptrAsyncParams->timesEta, ptrAsyncParams->vecEta);

}

void DistributionHandlerImpl::setOmega(size_t aNState, const stdVectorXd &aTimes, const std::vector< eventMap_t > &aOmegas) {

	ptrAsyncParams->nState = aNState;
	bool areTimesDifferent = copyVectorWithCheck(aTimes, ptrAsyncParams->timesOmega);
	bool areRatesDifferent = copyVectorOfSparseTensorWithCheck(aNState, aOmegas, ptrAsyncParams->vecOmega);

	if(areTimesDifferent) {
		dirtyAsyncRateShifts = true;
		updateSchedulerOperationTimeChange();
	}
	if(areRatesDifferent) {
		dirtyOmega = true;
		updateSchedulerOperationRateChange();
	}

	printParamTensorDebug("Omega", ptrAsyncParams->timesOmega, ptrAsyncParams->vecOmega);

}

void DistributionHandlerImpl::setMassSpeciationEvents(const stdVectorXd &massSpeciationTimes, const stdMatrixXd &massSpeciationProb) {

	bool areTimesDifferent = copyVectorWithCheck(massSpeciationTimes, ptrSyncParams->massSpeciationTimes);
	bool areProbsDifferent = copyMatrixWithCheck(massSpeciationProb, ptrSyncParams->massSpeciationProb);

	if(areTimesDifferent) {
		dirtySyncEventTimes = true;
		updateSchedulerOperationSyncTimeChange();
	}
	if(areProbsDifferent) {
		updateSchedulerOperationRateChange();
	}
	printParamVectorDebug("Mass speciation", ptrSyncParams->massSpeciationTimes, ptrSyncParams->massSpeciationProb);

}

void DistributionHandlerImpl::setMassExtinctionEvents(const stdVectorXd &massExtinctionTimes, const stdMatrixXd &massExtinctionProb) {

	bool areTimesDifferent = copyVectorWithCheck(massExtinctionTimes, ptrSyncParams->massExtinctionTimes);
	bool areProbsDifferent = copyMatrixWithCheck(massExtinctionProb, ptrSyncParams->massExtinctionProb);

	if(areTimesDifferent) {
		dirtySyncEventTimes = true;
		updateSchedulerOperationSyncTimeChange();
	}
	if(areProbsDifferent) {
		updateSchedulerOperationRateChange();
	}
	printParamVectorDebug("Mass extinction", ptrSyncParams->massExtinctionTimes, ptrSyncParams->massExtinctionProb);

}


void DistributionHandlerImpl::setMassExtinctionStateChangeProb(const std::vector< stdMatrixXd> &massExtinctionStateChangeProb) {

	bool areProbsDifferent = copyVectorOfProbMatrixWithCheck(massExtinctionStateChangeProb, ptrSyncParams->massExtinctionStateChangeProb);

	if(areProbsDifferent) {
		updateSchedulerOperationRateChange();
	}

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] Mass extinction state change prob" << std::endl;
		for(size_t i=0; i< ptrSyncParams->massExtinctionStateChangeProb.size(); ++i) {
			ss << ptrSyncParams->massExtinctionStateChangeProb[i] << std::endl;
		}
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setMassSamplingEvents(const stdVectorXd &massSamplingTimes, const stdMatrixXd &massSamplingProb) {

	bool areTimesDifferent = copyVectorWithCheck(massSamplingTimes, ptrSyncParams->massSamplingTimes);
	bool areProbsDifferent = copyMatrixWithCheck(massSamplingProb, ptrSyncParams->massSamplingProb);

	if(areTimesDifferent) {
		dirtySyncEventTimes = true;
		updateSchedulerOperationSyncTimeChange();
	}
	if(areProbsDifferent) {
		updateSchedulerOperationRateChange();
	}

	printParamVectorDebug("Mass sampling", ptrSyncParams->massSamplingTimes, ptrSyncParams->massSamplingProb);
}

void DistributionHandlerImpl::setMassDestrSamplingEvents(const stdVectorXd &massDestrSamplingTimes, const stdMatrixXd &massDestrSamplingProb) {

	bool areTimesDifferent = copyVectorWithCheck(massDestrSamplingTimes, ptrSyncParams->massDestrSamplingTimes);
	bool areProbsDifferent = copyMatrixWithCheck(massDestrSamplingProb, ptrSyncParams->massDestrSamplingProb);

	if(areTimesDifferent) {
		dirtySyncEventTimes = true;
		updateSchedulerOperationSyncTimeChange();
	}
	if(areProbsDifferent) {
		updateSchedulerOperationRateChange();
	}

	printParamVectorDebug("Mass destr sampling", ptrSyncParams->massDestrSamplingTimes, ptrSyncParams->massDestrSamplingProb);
}

void DistributionHandlerImpl::setSyncMonitors(const std::vector< double > &aSynchMonitoring) {

	bool areTimesDifferent = copyVectorWithCheck(aSynchMonitoring, synchMonitoring);
	if(areTimesDifferent) updateSchedulerOperationTimeChange();

}

double DistributionHandlerImpl::computeLogLikelihood() {

	assert(ptrTree);
	assert(ptrData);
	assert(ptrData->isReady());

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	// Create the tensors and sync event
	updateParameters();
	updateSyncEvents();

	// Init scheduler
	if(scheduler == NULL || schedulerOperation == RESET) {
		scheduler.reset( new Likelihood::Scheduler::BaseScheduler(ptrTree) );
		scheduler->setRateShiftEvents(ptrAsyncParams);
		scheduler->setSynchronousEvents(ptrSyncEvents);
		scheduler->setMonitoringProbes(synchMonitoring);
		scheduler->defineAndSetRescalingEvents();
		dirtyApproximator = true;
	} else if(schedulerOperation == UPDATE) {
		if(dirtyAsyncRateShifts) { // Refresh rate shifts
			scheduler->removeRateShiftEvents();
			scheduler->removeRescalingEvents();

			scheduler->setRateShiftEvents(ptrAsyncParams);
			scheduler->defineAndSetRescalingEvents();

			dirtyAsyncRateShifts = false;
		}
		if(dirtySyncEventTimes) { // Refresh sync events time
			scheduler->removeSynchronousEvents();
			scheduler->removeRescalingEvents();

			scheduler->setSynchronousEvents(ptrSyncEvents);
			scheduler->setMonitoringProbes(synchMonitoring);
			scheduler->defineAndSetRescalingEvents();
			dirtySyncEventTimes = false;
		}
		// TODO UPDATE MONITORING PROB
	}

	assert(scheduler);
	Likelihood::Conditions::conditionalProbability_t condType = Likelihood::Conditions::intToConditionalProbabilityType(condProbType);
	Likelihood::Integrator::integrationScheme_t      intType  = Likelihood::Integrator::intToIntegratorType(integrationScheme);

	if(approximator == NULL || dirtyApproximator) {
		assert((size_t)approxVersion >= 0 && (size_t)approxVersion < Likelihood::Approximator::APPROXIMATOR_NAMES.size());
		if(approxVersion == SEQUENTIAL_OPTIMIZED) {
			approximator = (Likelihood::Approximator::Factory::createSequentialTemplateCPU(intType, condType, ptrData, scheduler, ptrSyncEvents, ptrTensors));
		} else if(approxVersion == SEQUENTIAL_BRANCHWISE) {
			approximator = (Likelihood::Approximator::Factory::createSequentialBranchwiseCPU(intType, condType, ptrData, scheduler, ptrSyncEvents, ptrTensors));
	#if defined(_OPENMP)
		} else if(approxVersion == PARALLEL_OPTIMIZED) {
			Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);
			Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
			//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
			approximator = (Likelihood::Approximator::Factory::createParallelOpenMPCPU(intType, condType, ptrData, scheduler, ptrSyncEvents, ptrTensors));
		} else if(approxVersion == PARALLEL_BRANCHWISE) {
			Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);
			Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
			//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
			approximator = (Likelihood::Approximator::Factory::createParallelBranchwiseCPU(intType, condType, ptrData, scheduler, ptrSyncEvents, ptrTensors));
	#endif
		} else if(approxVersion == AUTO_TUNING) {
			approximator = (Likelihood::Approximator::Factory::createAutoTuningCPU(intType, condType, ptrData, scheduler, ptrSyncEvents, ptrTensors));
		} else {
			assert(approximator != NULL && "The approximator requested is not available (openmp approximators are only available if you compiled with openmp (-fopenmp)).");
		}
		dirtyApproximator = false;
	}

	approximator->setConditionalProbabilityCompatibilityMode(compatibilityMode);

	approximator->setQuasistationaryFrequencyMode(useQuasistationaryFrequency);


	if(applyTreeLikCorrection) {
		approximator->enableTreeLikelihoodCorrection();
	} else {
		approximator->disableTreeLikelihoodCorrection();
	}

	approximator->setDefaultDeltaT(initDeltaT);

	approximator->setPriorStateProbability(rootPrior);

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	double logLik = approximator->approximateLogLikelihood();

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")

	schedulerOperation = NONE;


	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "=============================================================" << std::endl;
		double debugLik = debugLikelihoodEvaluation();
		ss.precision(10);
		ss << "DEBUG LIK (full reset) : " << std::scientific << debugLik << std::endl;
		ss << "DEFAULT LIK (update) : " << std::scientific << debugLik << std::endl;
		ss << "Lik difference : " << std::scientific << debugLik - logLik << std::fixed << std::endl;
		ss << "=============================================================" << std::endl;
		debugChoseOutputStream(ss.str());
		assert(fabs(debugLik - logLik) < 5.0*Likelihood::Approximator::BaseApproximator::DEFAULT_REL_TOLERANCE);
	}

	return logLik;
}

mapHistories_t DistributionHandlerImpl::drawHistory() {
	return drawMultipleHistories(1).front();
}

mapHistories_t DistributionHandlerImpl::drawHistoryAndComputeRates(std::vector<double>& averageLambda, std::vector<double>& averageMu, std::vector<double>& averagePhi, std::vector<double>& averageDelta, std::vector<long>& numChanges) {

	Likelihood::Approximator::StochasticMapping::StochasticMappingApproxSharedPtr stochasticMappingApprox(createStochasticMappingApprox());

	stochasticMappingApprox->setConditionalProbabilityCompatibilityMode(compatibilityMode);
	stochasticMappingApprox->setQuasistationaryFrequencyMode(useQuasistationaryFrequency);

	if(applyTreeLikCorrection) {
		stochasticMappingApprox->enableTreeLikelihoodCorrection();
	} else {
		stochasticMappingApprox->disableTreeLikelihoodCorrection();
	}

	stochasticMappingApprox->setDefaultDeltaT(initDeltaT);

	stochasticMappingApprox->setPriorStateProbability(rootPrior);

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	stochasticMappingApprox->approximateLogLikelihood();

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")

	mapHistories_t theHistory = stochasticMappingApprox->drawHistoryAndComputeRates(averageLambda, averageMu, averagePhi, averageDelta, numChanges);

	schedulerOperation = NONE;

	return theHistory;

}

mapHistories_t DistributionHandlerImpl::drawAncestralStates() {
	return drawMultipleAncestralStates(1).front();
}

Likelihood::Approximator::StochasticMapping::StochasticMappingApproxSharedPtr DistributionHandlerImpl::createStochasticMappingApprox() {

	updateParameters();
	updateSyncEvents();

	// Init scheduler
	if(scheduler == NULL || schedulerOperation == RESET) {
		scheduler.reset( new Likelihood::Scheduler::BaseScheduler(ptrTree) );
		scheduler->setRateShiftEvents(ptrAsyncParams);
		scheduler->setSynchronousEvents(ptrSyncEvents);
		scheduler->setMonitoringProbes(synchMonitoring);
		scheduler->defineAndSetRescalingEvents();
		dirtyApproximator = true;
	} else if(schedulerOperation == UPDATE) {
		if(dirtyAsyncRateShifts) { // Refresh rate shifts
			scheduler->removeRateShiftEvents();
			scheduler->removeRescalingEvents();

			scheduler->setRateShiftEvents(ptrAsyncParams);
			scheduler->defineAndSetRescalingEvents();

			dirtyAsyncRateShifts = false;
		}
		if(dirtySyncEventTimes) { // Refresh sync events time
			scheduler->removeSynchronousEvents();
			scheduler->removeRescalingEvents();

			scheduler->setSynchronousEvents(ptrSyncEvents);
			scheduler->setMonitoringProbes(synchMonitoring);
			scheduler->defineAndSetRescalingEvents();
			dirtySyncEventTimes = false;
		}
		// TODO UPDATE MONITORING PROB
	}

	assert(scheduler);
	Likelihood::Conditions::conditionalProbability_t condType = Likelihood::Conditions::intToConditionalProbabilityType(condProbType);
	Likelihood::Integrator::integrationScheme_t      intType  = Likelihood::Integrator::intToIntegratorType(integrationScheme);

	return Likelihood::Approximator::Factory::createStochasticMappingApprox(intType, condType, ptrData, scheduler, ptrSyncEvents, ptrTensors);
}

vecHistories_t DistributionHandlerImpl::drawMultipleHistories(size_t nReplicas) {

	assert(nReplicas>0);

	vecHistories_t vecHistories;

	Likelihood::Approximator::StochasticMapping::StochasticMappingApproxSharedPtr stochasticMappingApprox(createStochasticMappingApprox());

	stochasticMappingApprox->setConditionalProbabilityCompatibilityMode(compatibilityMode);

	stochasticMappingApprox->setQuasistationaryFrequencyMode(useQuasistationaryFrequency);

	if(applyTreeLikCorrection) {
		stochasticMappingApprox->enableTreeLikelihoodCorrection();
	} else {
		stochasticMappingApprox->disableTreeLikelihoodCorrection();
	}

	stochasticMappingApprox->setDefaultDeltaT(initDeltaT);

	stochasticMappingApprox->setPriorStateProbability(rootPrior);

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	stochasticMappingApprox->approximateLogLikelihood();

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	for(size_t iR=0; iR<nReplicas; ++iR) {
		vecHistories.push_back(stochasticMappingApprox->drawHistory());
	}

	schedulerOperation = NONE;

	return vecHistories;
}

vecHistories_t DistributionHandlerImpl::drawMultipleAncestralStates(size_t nReplicas) {

	assert(nReplicas>0);

	vecHistories_t vecHistories;

	Likelihood::Approximator::StochasticMapping::StochasticMappingApproxSharedPtr stochasticMappingApprox(createStochasticMappingApprox());

	stochasticMappingApprox->setConditionalProbabilityCompatibilityMode(compatibilityMode);

	stochasticMappingApprox->setQuasistationaryFrequencyMode(useQuasistationaryFrequency);

	if(applyTreeLikCorrection) {
		stochasticMappingApprox->enableTreeLikelihoodCorrection();
	} else {
		stochasticMappingApprox->disableTreeLikelihoodCorrection();
	}

	stochasticMappingApprox->setDefaultDeltaT(initDeltaT);

	stochasticMappingApprox->setPriorStateProbability(rootPrior);

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	stochasticMappingApprox->approximateLogLikelihood();

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	for(size_t iR=0; iR<nReplicas; ++iR) {
		vecHistories.push_back(stochasticMappingApprox->drawAncestralStates());
	}

	schedulerOperation = NONE;

	return vecHistories;
}


size_t DistributionHandlerImpl::getVersion() const {
	return 1;
}

void DistributionHandlerImpl::setSeed(size_t aSeed) const {
	Utils::RNG::Manager::getInstance()->initializeDefaultRNG(aSeed);
}

void DistributionHandlerImpl::updateSchedulerOperationRateChange() {
	if(schedulerOperation == RESET) return;
	schedulerOperation = UPDATE;
}

void DistributionHandlerImpl::updateSchedulerOperationSyncTimeChange() {
	if(schedulerOperation == RESET) return;
	schedulerOperation = UPDATE;
}

void DistributionHandlerImpl::updateSchedulerOperationTimeChange() {
	schedulerOperation = RESET;
}


void DistributionHandlerImpl::updateParameters() {
	if(ptrTensors == NULL) {
		ptrTensors = Tensor::Factory::createContainerWithTimeHeterogenousVectors(ptrAsyncParams);
		dirtyLambda = dirtyMu = dirtyDelta = dirtyPhi = dirtyEta = dirtyOmega = false;
	} else {
		if(dirtyLambda) {
			ptrTensors->getLambda()->update();
			dirtyLambda = false;
		}
		if(dirtyMu) {
			ptrTensors->getMu()->update();
			dirtyMu = false;
		}
		if(dirtyDelta) {
			ptrTensors->getDelta()->update();
			dirtyDelta = false;
		}
		if(dirtyPhi) {
			ptrTensors->getPhi()->update();
			dirtyPhi = false;
		}
		if(dirtyDelta) {
			ptrTensors->getDelta()->update();
			dirtyDelta = false;
		}
		if(dirtyEta) {
			ptrTensors->getEta()->update();
			dirtyEta = false;
		}
		if(dirtyOmega) {
			ptrTensors->getOmega()->update();
			dirtyOmega = false;
		}
	}
}

void DistributionHandlerImpl::updateSyncEvents() {
	if(ptrSyncEvents == NULL) {
		ptrSyncEvents.reset(new SynchronousEvents::Container(ptrSyncParams));
	}
}


bool DistributionHandlerImpl::copyVectorWithCheck(const std::vector<double> &inVector, std::vector<double> &outVector) {

	bool hasChanged = false;
	if(outVector.size() != inVector.size()) {
		outVector.resize(inVector.size());
		hasChanged = true;
	}

	for(size_t i=0; i<inVector.size(); ++i) {
		hasChanged = hasChanged || (inVector[i] != outVector[i]);
		outVector[i] = inVector[i];
	}

	return hasChanged;
}

bool DistributionHandlerImpl::copyMatrixWithCheck(const stdMatrixXd &inMatrix, std::vector< Eigen::VectorXd > &outMatrix) {

	bool hasChanged = false;
	if(inMatrix.size() != outMatrix.size()) {
		outMatrix.resize(inMatrix.size());
		hasChanged = true;
	}

	for(size_t i=0; i<inMatrix.size(); ++i) {
		if(inMatrix[i].size() != (size_t)outMatrix[i].size()) {
			outMatrix[i].resize(inMatrix[i].size());
			hasChanged = true;
		}
		for(size_t j=0; j<inMatrix[i].size(); ++j) {
			hasChanged = hasChanged || (inMatrix[i][j] != outMatrix[i][j]);
			outMatrix[i](j) = inMatrix[i][j];
		}
	}

	return hasChanged;
}


bool DistributionHandlerImpl::copyVectorOfRateMatrixWithCheck(const std::vector< stdMatrixXd> &inVecMatrix, std::vector< Eigen::MatrixXd > &outVecMatrix) {

	bool hasChanged = false;

	if(inVecMatrix.size() != outVecMatrix.size()) {
		outVecMatrix.resize(inVecMatrix.size());
		hasChanged = true;
	}

	for(size_t iV=0; iV<inVecMatrix.size(); ++iV) {
		// Checking size change
		if(inVecMatrix[iV].size() != (size_t)outVecMatrix[iV].cols() || inVecMatrix[iV].size() != (size_t)outVecMatrix[iV].rows()) {
			outVecMatrix[iV].resize(inVecMatrix[iV].size(), inVecMatrix[iV].size());
			hasChanged = true;
		}

		// Copying and checking for rate change
		for(size_t iX=0; iX<inVecMatrix[iV].size(); ++iX) {
			double rowSum = 0;
			for(size_t iY=0; iY<inVecMatrix[iV][iX].size(); ++iY) {
				hasChanged = hasChanged || (inVecMatrix[iV][iX][iY] != outVecMatrix[iV](iX, iY));
				outVecMatrix[iV](iX, iY) = inVecMatrix[iV][iX][iY];
				rowSum += inVecMatrix[iV][iX][iY];
			}

			// Checking if this matrix make sense
			if(rowSum != 0. && fabs(rowSum) < 1e-7) {
				outVecMatrix[iV](iX, iX) = 0.;
				outVecMatrix[iV](iX, iX) -= outVecMatrix[iV].row(iX).sum();
			} else if (rowSum != 0.) {
				assert(fabs(rowSum) >= 1e-7 && "A row from an Eta matrix does not sum to 0.");
			}

		}
	}

	return hasChanged;
}

bool DistributionHandlerImpl::copyVectorOfProbMatrixWithCheck(const std::vector< stdMatrixXd>  &inVecMatrix, std::vector< Eigen::MatrixXd > &outVecMatrix) {

	bool hasChanged = false;

	if(inVecMatrix.size() != outVecMatrix.size()) {
		outVecMatrix.resize(inVecMatrix.size());
		hasChanged = true;
	}

	for(size_t iV=0; iV<inVecMatrix.size(); ++iV) {
		// Checking size change
		if(inVecMatrix[iV].size() != (size_t)outVecMatrix[iV].cols() || inVecMatrix[iV].size() != (size_t)outVecMatrix[iV].rows()) {
			outVecMatrix[iV].resize(inVecMatrix[iV].size(), inVecMatrix[iV].size());
			hasChanged = true;
		}

		// Copying and checking for rate change
		for(size_t iX=0; iX<inVecMatrix[iV].size(); ++iX) {
			double rowSum = 0;
			for(size_t iY=0; iY<inVecMatrix[iV].size(); ++iY) {
				hasChanged = hasChanged || (inVecMatrix[iV][iX][iY] != outVecMatrix[iV](iX, iY));
				outVecMatrix[iV](iX, iY) = inVecMatrix[iV][iX][iY];
				rowSum += inVecMatrix[iV][iX][iY];
			}

			// Checking if this matrix make sense
			assert(fabs(1.-rowSum) < 1e-10 && "A row from a probability matrix does not sum to 0.");

		}
	}

	return hasChanged;
}


bool DistributionHandlerImpl::copyVectorOfSparseTensorWithCheck(size_t aN, const std::vector< eventMap_t> &inVecSpTensors, std::vector< eventMap_t> &outVecSpTensors) {

	bool hasChanged = true;

	outVecSpTensors = inVecSpTensors;

	for(size_t iV=0; iV<outVecSpTensors.size(); ++iV) {

		// Copy and sum over the dimensions
		std::vector<double> sums(aN, 0.);
		for(eventMap_t::const_iterator it = outVecSpTensors[iV].begin(); it != outVecSpTensors[iV].end(); ++it) {

			std::vector<unsigned int> key = it->first;
			double val = it->second;

			assert(key.size() == 3 && key[0] < aN && key[1] < aN && key[2] < aN);

			sums[key[0]] += val;
		}

		for(size_t iS=0; iS < sums.size(); ++iS) {
			// std::cout << "row " << iS << " sum to " << sums[iS] << std::endl;
			assert(fabs(1. - sums[iS]) < 1.e-10 && "Each row should sum to one.");
		}
	}

	return hasChanged;
}


void DistributionHandlerImpl::printParamVectorDebug(const std::string &paramName, const std::vector<double> &aTime, const std::vector< Eigen::VectorXd > &params) {
	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] " << paramName << std::endl;
		ss << "Times = " ;
		for(size_t iT=0; iT<aTime.size(); ++iT) ss << aTime[iT] << " , ";
		ss << std::endl;
		for(size_t iP = 0; iP < params.size(); ++iP) {
			ss << "Param [ " << iP << " ] " << params[iP].transpose() << std::endl;
		}
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}
}

void DistributionHandlerImpl::printParamMatrixDebug(const std::string &paramName, const std::vector<double> &aTime, const std::vector< Eigen::MatrixXd > &params) {
	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] " << paramName << std::endl;
		ss << "Times = " ;
		for(size_t iT=0; iT<aTime.size(); ++iT) ss << aTime[iT] << " , ";
		ss << std::endl;
		for(size_t iP = 0; iP < params.size(); ++iP) {
			ss << "Param [ " << iP << " ] " << std::endl << params[iP] << std::endl;
		}
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}
}

void DistributionHandlerImpl::printParamTensorDebug(const std::string &paramName, const std::vector<double> &aTime, const std::vector< eventMap_t > &params) {

	if(debugMode != DBG_NONE) {
		std::stringstream ss;
		ss << "[DEBUG TENSORPHYLO] " << paramName << std::endl;
		ss << "Times = " ;
		for(size_t iT=0; iT<aTime.size(); ++iT) ss << aTime[iT] << " , ";
		ss << std::endl;
		for(size_t iP = 0; iP < params.size(); ++iP) {
			ss << "Param [ " << iP << " ] " << std::endl ;
			for(eventMap_t::const_iterator it=params[iP].begin(); it != params[iP].end(); ++it) {
				ss << it->first[0] << ", " << it->first[1] << ", " << it->first[2] << " : " << it->second << std::endl;
			}
		}
		ss << "-----------------------------------------------------" << std::endl;
		debugChoseOutputStream(ss.str());
	}

}

void DistributionHandlerImpl::setDebugMode(debugMode_t aDebugMode) {
	assert(aDebugMode != DBG_FILE);
	debugMode = aDebugMode;
}

void DistributionHandlerImpl::setDebugMode(debugMode_t aDebugMode, const std::string &aFilePath) {
	assert(aDebugMode == DBG_FILE);

	debugMode = aDebugMode;
	debugFilePath = aFilePath;

	if(debugFile.is_open()) debugFile.close();

	debugFile.open(aFilePath.c_str());
	assert(debugFile.good());
}

void DistributionHandlerImpl::writeStateToFile(const std::string &aFilePath) {

	std::ofstream myFile;
	myFile.open(aFilePath.c_str(), std::ofstream::out);
	assert(myFile.good());

	myFile << "applyTreeLikCorrection" << std::endl;
	if(applyTreeLikCorrection) {
		myFile << "true" << std::endl << std::endl;
	} else {
		myFile << "false" << std::endl << std::endl;
	}

	myFile << "useQuasistationaryFrequency" << std::endl;
	if(useQuasistationaryFrequency) {
		myFile << "true" << std::endl << std::endl;
	} else {
		myFile << "false" << std::endl << std::endl;
	}

	myFile << "compatibilityMode" << std::endl;
	if(compatibilityMode) {
		myFile << "true" << std::endl << std::endl;
	} else {
		myFile << "false" << std::endl << std::endl;
	}

	myFile << "intLikApproximator" << std::endl;
	myFile << static_cast<int>(approxVersion) << std::endl << std::endl;

	myFile << "nThreads" << std::endl;
	myFile << nThreads << std::endl << std::endl;

	myFile << "intCondType" << std::endl;
	myFile << static_cast<int>(condProbType) << std::endl << std::endl;

	myFile << "deltaT" << std::endl;
	myFile << initDeltaT << std::endl << std::endl;

	myFile << "rootPrior" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(rootPrior) << std::endl ;

	myFile << "timesLambda" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->timesLambda) << std::endl;

	myFile << "lambda" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->vecLambda) << std::endl;

	myFile << "timesMu" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->timesMu) << std::endl;

	myFile << "mu" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->vecMu) << std::endl;

	myFile << "timesPhi" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->timesPhi) << std::endl;

	myFile << "phi" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->vecPhi) << std::endl;

	myFile << "timesDelta" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->timesDelta) << std::endl;

	myFile << "delta" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->vecDelta) << std::endl;

	myFile << "timesEta" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->timesEta) << std::endl;

	myFile << "eta" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->vecEta) << std::endl;

	myFile << "timesOmega" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->timesOmega) << std::endl;

	myFile << "omega" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrAsyncParams->vecOmega) << std::endl;

	myFile << "massSpeciationTimes" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massSpeciationTimes) << std::endl;

	myFile << "massSpeciationProb" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massSpeciationProb) << std::endl;

	myFile << "massExtinctionTimes" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massExtinctionTimes) << std::endl;

	myFile << "massExtinctionProb" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massExtinctionProb) << std::endl;

	myFile << "massExtinctionStateChangeProb" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massExtinctionStateChangeProb) << std::endl;

	myFile << "massSamplingTimes" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massSamplingTimes) << std::endl;

	myFile << "massSamplingProb" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massSamplingProb) << std::endl;

	myFile << "massDestrSamplingTimes" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massDestrSamplingTimes) << std::endl;

	myFile << "massDestrSamplingProb" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(ptrSyncParams->massDestrSamplingProb) << std::endl;

	myFile << "synchMonitoring" << std::endl;
	myFile << Utils::Serializer::Parameters::toString(synchMonitoring) << std::endl;

	myFile.close();

}

void DistributionHandlerImpl::loadStateFromFile(const std::string &aFilePath) {

	std::ifstream myFile;
	myFile.open(aFilePath.c_str(), std::ifstream::in);
	assert(myFile.good());

	std::string line;
	while(std::getline(myFile, line)) {

		boost::algorithm::trim(line);
		if(line.empty()) continue;

		if(line =="applyTreeLikCorrection") {
			std::string strApplyTLC;
			myFile >> strApplyTLC;
			boost::algorithm::to_lower<std::string>(strApplyTLC);
			assert((strApplyTLC == "true" || strApplyTLC == "false") && "The value for 'likTreeCorrection' must be 'true' or 'false'.");
			applyTreeLikCorrection = strApplyTLC == "true";
		} else if(line =="useQuasistationaryFrequency") {
			std::string strQuasiFrequencyMode;
			myFile >> strQuasiFrequencyMode;
			boost::algorithm::to_lower<std::string>(strQuasiFrequencyMode);
			assert((strQuasiFrequencyMode == "true" || strQuasiFrequencyMode == "false") && "The value for 'QuasistationaryFrequency Mode' must be 'true' or 'false'.");
			useQuasistationaryFrequency = strQuasiFrequencyMode == "true";
		} else if(line =="compatibilityMode") {
			std::string strCompaMode;
			myFile >> strCompaMode;
			boost::algorithm::to_lower<std::string>(strCompaMode);
			assert((strCompaMode == "true" || strCompaMode == "false") && "The value for 'Compatibility Mode' must be 'true' or 'false'.");
			compatibilityMode = strCompaMode == "true";
		} else if(line =="intLikApproximator") {
			int intLikApproximator = 0;
			myFile >> intLikApproximator;
			approxVersion = static_cast<approximatorVersion_t>(intLikApproximator);
		} else if(line == "nThreads") {
			myFile >> nThreads;
		} else if(line =="intCondType") {
			int intCondType = 0;
			myFile >> intCondType;
			condProbType = static_cast<conditionalProbability_t>(intCondType);;
		} else if(line =="deltaT") {
			myFile >> initDeltaT;
		} else if(line =="rootPrior") {
			Utils::Serializer::Parameters::fromFile(myFile, rootPrior);
		} else if(line =="timesLambda") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->timesLambda);
		} else if(line =="lambda") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->vecLambda);
		} else if(line =="timesMu") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->timesMu);
		} else if(line =="mu") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->vecMu);
		} else if(line =="timesPhi") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->timesPhi);
		} else if(line =="phi") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->vecPhi);
		} else if(line =="timesDelta") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->timesDelta);
		} else if(line =="delta") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->vecDelta);
		} else if(line =="timesEta") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->timesEta);
		} else if(line =="eta") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->vecEta);
		} else if(line =="timesOmega") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->timesOmega);
		} else if(line =="omega") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrAsyncParams->vecOmega);
		} else if(line =="massSpeciationTimes") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massSpeciationTimes);
		} else if(line =="massSpeciationProb") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massSpeciationProb);
		} else if(line =="massExtinctionTimes") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massExtinctionTimes);
		} else if(line =="massExtinctionProb") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massExtinctionProb);
		} else if(line =="massExtinctionStateChangeProb") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massExtinctionStateChangeProb);
		} else if(line =="massSamplingTimes") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massSamplingTimes);
		} else if(line =="massSamplingProb") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massSamplingProb);
		} else if(line =="massDestrSamplingTimes") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massDestrSamplingTimes);
		} else if(line =="massDestrSamplingProb") {
			Utils::Serializer::Parameters::fromFile(myFile, ptrSyncParams->massDestrSamplingProb);
		} else if(line =="synchMonitoring") {
			Utils::Serializer::Parameters::fromFile(myFile, synchMonitoring);
		} else {
			std::cerr << "Parameters reader error at string: " << line << std::endl;
		}

	}

	myFile.close();

	// force reinit:
	Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
	Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);

	// Create the tensors and sync event
	ptrTensors.reset();
	updateParameters();
	ptrSyncEvents.reset();
	updateSyncEvents();

	schedulerOperation = RESET;
	dirtyAsyncRateShifts = dirtySyncEventTimes = dirtyApproximator = false;
	dirtyLambda = dirtyMu = dirtyDelta = dirtyPhi = dirtyEta = dirtyOmega = false;
}

double DistributionHandlerImpl::debugLikelihoodEvaluation() {

	assert(ptrTree);
	assert(ptrData);
	assert(ptrData->isReady());

	//Utils::Output::outputManager().setVerbosityThreshold(Utils::Output::MEDIUM_VERB);

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
	//Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);

	// Create the tensors and sync event
	Tensor::ContainerSharedPtr ptrDbgTensors(Tensor::Factory::createContainerWithTimeHeterogenousVectors(ptrAsyncParams));
	SynchronousEvents::ContainerSharedPtr ptrDbgSyncEvents(new SynchronousEvents::Container(ptrSyncParams));

	// Init scheduler
	//ss << ptrTree->getNewickString() << std::endl;
	boost::shared_ptr<Likelihood::Scheduler::BaseScheduler> ptrDbgScheduler( new Likelihood::Scheduler::BaseScheduler(ptrTree) );

	ptrDbgScheduler->setRateShiftEvents(ptrAsyncParams);
	ptrDbgScheduler->setSynchronousEvents(ptrDbgSyncEvents);
	ptrDbgScheduler->setMonitoringProbes(synchMonitoring);
	ptrDbgScheduler->defineAndSetRescalingEvents();

	Likelihood::Conditions::conditionalProbability_t condType = Likelihood::Conditions::intToConditionalProbabilityType(condProbType);
	Likelihood::Integrator::integrationScheme_t      intType  = Likelihood::Integrator::intToIntegratorType(integrationScheme);
	Likelihood::Approximator::ApproximatorSharedPtr ptrDbgApprox;

	assert(ptrDbgScheduler);

	assert((size_t)approxVersion >= 0 && (size_t)approxVersion < Likelihood::Approximator::APPROXIMATOR_NAMES.size());
	if(approxVersion == SEQUENTIAL_OPTIMIZED) {
		ptrDbgApprox = (Likelihood::Approximator::Factory::createAutoTuningCPU(intType, condType, ptrData, ptrDbgScheduler, ptrDbgSyncEvents, ptrDbgTensors));
	} else if(approxVersion == SEQUENTIAL_OPTIMIZED) {
		ptrDbgApprox = (Likelihood::Approximator::Factory::createSequentialTemplateCPU(intType, condType, ptrData, ptrDbgScheduler, ptrDbgSyncEvents, ptrDbgTensors));
	} else if(approxVersion == SEQUENTIAL_BRANCHWISE) {
		ptrDbgApprox = (Likelihood::Approximator::Factory::createSequentialBranchwiseCPU(intType, condType, ptrData, ptrDbgScheduler, ptrDbgSyncEvents, ptrDbgTensors));
#if defined(_OPENMP)
	} else if(approxVersion == PARALLEL_OPTIMIZED) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
		//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
		ptrDbgApprox = (Likelihood::Approximator::Factory::createParallelOpenMPCPU(intType, condType, ptrData, ptrDbgScheduler, ptrDbgSyncEvents, ptrDbgTensors));
	} else if(approxVersion == PARALLEL_BRANCHWISE) {
		Utils::Parallel::Manager::getInstance()->setMaxNThread(nThreads);
		Utils::Parallel::Manager::getInstance()->setNThread(nThreads);
		//assert(Utils::Parallel::Manager::getInstance()->useOpenMP() && "This code is not compiled with openmp (-fopenmp) or you requested to use only 1 processor (use approximator 'Likelihood::Approximator::SEQUENTIAL_OPTIMIZED'.");
		ptrDbgApprox = (Likelihood::Approximator::Factory::createParallelBranchwiseCPU(intType, condType, ptrData, ptrDbgScheduler, ptrDbgSyncEvents, ptrDbgTensors));
#endif
	} else if(approxVersion == AUTO_TUNING) {
		ptrDbgApprox = (Likelihood::Approximator::Factory::createAutoTuningCPU(intType, condType, ptrData, ptrDbgScheduler, ptrDbgSyncEvents, ptrDbgTensors));
	} else {
		assert(ptrDbgApprox != NULL && "The approximator requested is not available (openmp approximators are only available if you compiled with openmp (-fopenmp)).");
	}

	assert(ptrDbgApprox);

	ptrDbgApprox->setConditionalProbabilityCompatibilityMode(compatibilityMode);

	ptrDbgApprox->setQuasistationaryFrequencyMode(useQuasistationaryFrequency);

	if(applyTreeLikCorrection) {
		ptrDbgApprox->enableTreeLikelihoodCorrection();
	} else {
		ptrDbgApprox->disableTreeLikelihoodCorrection();
	}

	ptrDbgApprox->setDefaultDeltaT(initDeltaT);

	ptrDbgApprox->setPriorStateProbability(rootPrior);

	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "INIT_TIME")

	//_START_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")
	double logLik = ptrDbgApprox->approximateLogLikelihood();
	//_END_EVENT(Utils::Profiling::UniqueProfiler::getInstance(), "TIME_AFTER_INIT")

	return logLik;
}

void DistributionHandlerImpl::debugChoseOutputStream(const std::string &string) {
	if(debugMode == DBG_PRINT) {
		std::cout << string;
		std::cout.flush();
	} else if(debugMode == DBG_FILE) {
		debugFile << string;
		debugFile.flush();
	} else {
		assert(false && "Unknown output stream.");
	}
}

Eigen::VectorXd DistributionHandlerImpl::getQuasiStationaryFrequency(double t) {

	// update parameters
	updateParameters();

	// just call the internal function
	return Likelihood::Kernels::CPU::computeQuasiStationaryFrequency(ptrTensors, t);

}


} /* namespace Interface */
} /* namespace TensorPhylo */
