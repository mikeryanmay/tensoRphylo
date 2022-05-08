/*
 * SequentialStochasticMappingCPU.cpp
 *
 *  Created on: May 6, 2020
 *      Author: xaviermeyer
 */

#include "SequentialStochasticMappingCPU.h"

#include <cmath>
#include <algorithm>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>


#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "../IncLikelihoodApproximator.h"
#include "../Specialized/AdaptiveSegmentFW.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Utils/Output/OutputManager.h"
#include "Likelihood/Scheduler/DAG/IncDAG.h"
#include "Likelihood/CustomIntegrators/AdaptiveIntegrators.hpp"
#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "Likelihood/StateTypes/Vector/EigenStateOperations.hpp"
#include "Likelihood/Kernels/CPU/Specialized/EigenKernelsSingleton.h"
#include "Likelihood/Kernels/CPU/Specialized/IntegrationKernelSingleton.h"
#include "../Specialized/DenseSegmentBW.h"
#include "../Specialized/SingleStepSegmentFW.h"
#include "SegmentStatistics.h"
#include "Utils/RNG/Manager.h"
#include "Utils/RNG/Sitmo11RNG.h"

using Utils::RNG::Sitmo11RNG;

namespace Likelihood {
namespace Approximator {
namespace StochasticMapping {

const double SequentialStochasticMappingCPU::N_STEPS_PER_TREE = 500.;
const size_t SequentialStochasticMappingCPU::REJECTION_SAMPLING_TRY_LIMIT = 1e6;
const double SequentialStochasticMappingCPU::PROB_THRESHOLD = 0.02;
const double SequentialStochasticMappingCPU::PROB_TOLERANCE = PROB_THRESHOLD*1e-1;

SequentialStochasticMappingCPU::SequentialStochasticMappingCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
	   	   	   	   	   	   	    Conditions::conditionalProbability_t aConditionType,
								Phylogeny::Data::ContainerSharedPtr aPtrData,
								Scheduler::SchedulerSharedPtr aPtrScheduler,
								SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
								Tensor::ContainerSharedPtr aPtrTensorCont) :
		AsynchronousApproximator(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
		BaseStochasticMappingApproximator(),
		FIXED_EULER_DT(ptrScheduler->getPtrTree()->getOldestNode()->getAge()/N_STEPS_PER_TREE),
		kernelsBW(aConditionType, ptrData, ptrSApproxManager, ptrSyncEventsCont, ptrTensorCont, condCompatibilityMode),
		intKernelBW(ptrSApproxManager, ptrTensorCont),
		intKernelFW(ptrSApproxManager, ptrTensorCont){

	accAdaptiveDT(-0.1);
	stochasticMappingAlgo = REJECTION_SAMPLING_ALGO;
//	stochasticMappingAlgo = DENSE_EULER_ALGO;
//	stochasticMappingAlgo = DENSE_DOPRI_ALGO;

	// Important: initialize the size of vector that the memory pull will allocate.
	Utils::MemoryPool::eigenCPU().setNCategories(ptrTensorCont->getNumberOfState());

	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<stateType_t, intKernelBW_t, operations_t >(DEFAULT_ABS_TOLERANCE, DEFAULT_REL_TOLERANCE, deltaT, aIntScheme);
}

SequentialStochasticMappingCPU::~SequentialStochasticMappingCPU() {
	delete ptrIntegrator;
}

void SequentialStochasticMappingCPU::setAlgorithm(stochasticMappingAlgo_t aAlgo) {
	stochasticMappingAlgo = aAlgo;
}

void SequentialStochasticMappingCPU::setDefaultDeltaT(double aDeltaT) {
	assert(ptrIntegrator != NULL);
	deltaT = aDeltaT;
	ptrIntegrator->setDeltaT(deltaT);
}

size_t SequentialStochasticMappingCPU::getTotalNumberOfIntegrationSteps() const {
	return ptrIntegrator->getNSteps();
}

void SequentialStochasticMappingCPU::doPreProcessingSteps() {
	// Init with an empty state
	probState = Likelihood::StateType::Vector::EigenState();

	if(ptrScheduler->hasBeenUpdated()) {
		// reset DAG
		ptrDAG->rebuild();
		ptrScheduler->clearHasBeenUpdatedFlag();
	}

	ptrSchedDAG->reset();

	integrationTimes.clear();

	checkpoints.clear();
	tasksStack.clear();

	vecProbesState.clear();

}

void SequentialStochasticMappingCPU::doIntegrationStep(Likelihood::Scheduler::DAG::NodeDAG* task) {
	// We can get the layer of edges using the scheduler
	// However, it's not even needed as the kernels already deal with that
	// After a doEventStep, the state has already been resized to fit the next set of edges

	mapStatesVal_t statesVal;

	// iEdgesLayer starts at iEvent = iEdgesLayer and ends at iEvent+1=iEdgesLayer+1
	assert(task->isIntegrationTask());
	double startTime = task->getStartTime();
	double endTime = task->getEndTime();



	if(endTime == startTime) {
		statesVal.first = *task->getProbStates().front();
		statesVal.second = *task->getProbStates().front();
		checkpoints[task->getId()] = statesVal;
		return;
	}

	// We make sure that endTime is always within -1*epsilon from the true value
	// to make sure that we can find in which interval we are
	endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	if(endTime != ptrScheduler->getEvents().back()->getTime()) {
		endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	}

	assert(task->getProbStates().size() == 1);
	stateType_t &activePState = *task->getProbStates().front();
	statesVal.first = activePState;

	assert(ptrIntegrator != NULL);
	ptrIntegrator->integrate(startTime, endTime, activePState, intKernelBW);

	integrationTimes.insert(integrationTimes.end(), ptrIntegrator->getVecTimes().begin(), ptrIntegrator->getVecTimes().end());

	statesVal.second = activePState;
	checkpoints[task->getId()] = statesVal;

	tasksStack.push_back(task);
}

void SequentialStochasticMappingCPU::doEventStep(Likelihood::Scheduler::DAG::NodeDAG* task) {

	Likelihood::Scheduler::Event* event = task->getEvent();

	if(event->checkEvent(Likelihood::Scheduler::PRESENT_TIME_EVENT)) {			// When we get a sampled ancestor, we need to be careful when we are in the single branch hack case
		kernelsBW.setInitialCondition(task);
	} else if(event->checkEvent(Likelihood::Scheduler::NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isSpeciationNode()) {
			kernelsBW.computeAsynchSpeciation(event->getTime(), task);
		} else if(eventNode->isSampledAncestor()) {
			// When we get a sampled ancestor, we need to be careful when we are in the single branch hack case
			kernelsBW.computeAsynchSampling(event->getTime(), task);
		} else if(eventNode->isExtinct()) {
			kernelsBW.setInitalExtinctNodeCondition(event->getTime(), task);
		}
	} else if(event->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isRootNode() && eventNode->isSpeciationNode()) {
			kernelsBW.computeAsynchSpeciation(event->getTime(), task);
		}
		Eigen::VectorXd rf = this->getRootFrequency(event->getTime());
		logLikelihood = kernelsBW.computeLogLikelihood(event->getTime(), rf, task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT)) {
		kernelsBW.computeMassSpeciation(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT)) {
		kernelsBW.computeMassExtinction(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT)) {
		kernelsBW.computeMassSamplingEvent(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		kernelsBW.computeMassDestrSamplingEvent(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RATE_SHIFT)) {
		// DO NOTHING BUT ENSURE THAT THE NUMBERICAL INTEGRATOR DEAL ONLY WITH CONTIUNOUS FUNCTIONS
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RESCALING_EVENT)) {
		kernelsBW.rescaleProbabilities(task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		kernelsBW.rescaleProbabilities(task); // Forcing rescaling before probe
		this->doReportState(task);
	} else {
		assert(false && "Event is not implemented.");
	}

	if(Utils::Output::outputManager().check(Utils::Output::HIGH_VERB) && !event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		std::cout << "Event :" << event->toString() << std::endl;
		std::cout << "State : " << probState.toString();
		std::cout << "----------------------------------------" << std::endl;
	}

	tasksStack.push_back(task);

}

void SequentialStochasticMappingCPU::doPostProcessingSteps() {
	// Nothing to be done?
	assert(ptrDAG->getRoot()->getProbStates().size() == 1);
	stateType_t *ptrRootState = ptrDAG->getRoot()->getProbStates().front();
	probState = *ptrRootState; // Copy

	// Cleanup
	ptrRootState->removeVecProb();
	ptrDAG->getRoot()->detachStates();
	delete ptrRootState;

	orderProbes();

	assert(ptrDAG->getRoot()->getProbStates().empty());
}

void SequentialStochasticMappingCPU::doReportState(Likelihood::Scheduler::DAG::NodeDAG* task) {

	if(task->getProbStates().empty()) return;

	Likelihood::Monitor::ProbeState probe;
	probe.time = task->getTime();
	ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(probe.time, probe.u);

	for(size_t iP=0; iP<task->getProbStates().size(); ++iP) {
		probe.vecP.push_back(task->getProbStates()[iP]->getStateProb());
		probe.vecIdEdge.push_back(task->getProbStates()[iP]->getEdgeMapping());
		probe.scalingFactor.push_back(task->getProbStates()[iP]->getScaling());
	}

	std::pair<double, bool> resSR = intKernelBW.getStiffnessRatio(probe.time, probe.u);
	probe.stiffnessRatio = resSR.first;
	probe.hasNegativeImgEIGVal = resSR.second;

	vecProbesState.push_back(probe);
}

mapHistories_t SequentialStochasticMappingCPU::drawHistory() {

	mapHistories_t histories;

	for(long int iT = tasksStack.size()-1; iT >= 0; --iT) {
		Likelihood::Scheduler::DAG::NodeDAG* task = tasksStack[iT];

		if(task->isEventTask()) {
			doSimulateEvent(task, histories);
		} else if(task->isIntegrationTask()) {
			doSimulateBranch(task, histories);
		} else {
			assert(false && "Unknown task type");
		}
	}

	// Branch cleanup
	for(mapHistories_t::iterator itHistory = histories.begin(); itHistory != histories.end(); ++itHistory) {

		std::vector< std::pair<double, size_t> > branchHistory = itHistory->second;
		itHistory->second.clear();
		double lastTransitionAge = branchHistory.front().first;
		size_t lastTransitionState = branchHistory.front().second;
		for(size_t iT=1; iT<branchHistory.size(); ++iT) {
			if(branchHistory[iT].second != lastTransitionState || iT == branchHistory.size()-1) {
				double duration = lastTransitionAge - branchHistory[iT].first;
				itHistory->second.push_back(std::make_pair(duration, lastTransitionState));
				lastTransitionAge = branchHistory[iT].first;
				lastTransitionState = branchHistory[iT].second;
			}
		}
	}
	return histories;
}

mapHistories_t SequentialStochasticMappingCPU::drawHistoryAndComputeRates(std::vector<double>& averageLambda, std::vector<double>& averageMu, std::vector<double>& averagePhi, std::vector<double>& averageDelta, std::vector<long>& numChanges) {

	mapHistories_t histories;

	for(long int iT = tasksStack.size()-1; iT >= 0; --iT) {
		Likelihood::Scheduler::DAG::NodeDAG* task = tasksStack[iT];

		if(task->isEventTask()) {
			doSimulateEvent(task, histories);
		} else if(task->isIntegrationTask()) {
			doSimulateBranch(task, histories);
		} else {
			assert(false && "Unknown task type");
		}
	}

	// Make sure the branch-rate containers are clear
	size_t num_histories = histories.size() - 1;
	averageLambda.resize(num_histories);
	averageMu.resize(num_histories);
	averagePhi.resize(num_histories);
	averageDelta.resize(num_histories);
	numChanges.resize(num_histories);

	// Branch cleanup
	int avg_index = 0;
	for(mapHistories_t::iterator itHistory = histories.begin(); itHistory != histories.end(); ++itHistory) {

		if ( itHistory == histories.begin() ) {
			// skip the first history (it's the stem)
			continue;
		}

		std::vector< std::pair<double, size_t> > branchHistory = itHistory->second;
		itHistory->second.clear();
		double lastTransitionAge = branchHistory.front().first;
		size_t lastTransitionState = branchHistory.front().second;

		double totalBranchLength = 0.0;
		double totalLambda = 0.0;
		double totalMu = 0.0;
		double totalPhi = 0.0;
		double totalDelta = 0.0;
		long   totalNumChanges = 0;

		for(size_t iT=1; iT<branchHistory.size(); ++iT) {

//			if(branchHistory[iT].second != lastTransitionState || iT == branchHistory.size()-1) {

				// add the transition to the history
				double duration = lastTransitionAge - branchHistory[iT].first;
				itHistory->second.push_back(std::make_pair(duration, lastTransitionState));

				// increment the branch length
				totalBranchLength += duration;

				// increment the branch rates
				totalLambda += (ptrTensorCont->getEigenVecLambda(lastTransitionAge).col(0))(lastTransitionState) * duration;
				totalMu     += (ptrTensorCont->getEigenVecMu(lastTransitionAge).col(0))(lastTransitionState) * duration;
				totalPhi    += (ptrTensorCont->getEigenVecPhi(lastTransitionAge).col(0))(lastTransitionState) * duration;
				totalDelta  += (ptrTensorCont->getEigenVecDelta(lastTransitionAge).col(0))(lastTransitionState) * duration;

				if (branchHistory[iT].second != lastTransitionState)
				{
					totalNumChanges += 1;
				}

				// reset the event
				lastTransitionAge = branchHistory[iT].first;
				lastTransitionState = branchHistory[iT].second;

//			}

		} // end loop over history

		// compute the average rates
		averageLambda[avg_index] = totalLambda / totalBranchLength;
		averageMu[avg_index]     = totalMu / totalBranchLength;
		averagePhi[avg_index]    = totalPhi / totalBranchLength;
		averageDelta[avg_index]  = totalDelta / totalBranchLength;
		numChanges[avg_index++]  = totalNumChanges;

		// averageLambda[itHistory->first - 1] = totalLambda / totalBranchLength;
		// averageMu[itHistory->first - 1]     = totalMu / totalBranchLength;
		// averagePhi[itHistory->first - 1]    = totalPhi / totalBranchLength;
		// averageDelta[itHistory->first - 1]  = totalDelta / totalBranchLength;
		// numChanges[itHistory->first - 1]    = totalNumChanges;

	}

	return histories;

}

mapHistories_t SequentialStochasticMappingCPU::drawAncestralStates() {

	mapHistories_t histories;

	for(long int iT = tasksStack.size()-1; iT >= 0; --iT) {
		Likelihood::Scheduler::DAG::NodeDAG* task = tasksStack[iT];

		if(task->isEventTask()) {
			doSimulateEvent(task, histories);
		} else if(task->isIntegrationTask()) {
			doDrawAncestralState(task, histories);
		} else {
			assert(false && "Unknown task type");
		}
	}
	// Branch cleanup
	for(mapHistories_t::iterator itHistory = histories.begin(); itHistory != histories.end(); ++itHistory) {
		std::vector< std::pair<double, size_t> > branchHistory = itHistory->second;
		itHistory->second.clear();

		itHistory->second.push_back(branchHistory.front());
		itHistory->second.push_back(branchHistory.back());
	}

	return histories;
}

void SequentialStochasticMappingCPU::doSimulateAsyncSpeciationEvent(size_t parentState, Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {
	Likelihood::Scheduler::Event* event = task->getEvent();

	// Recover children subtree conditional prob
	assert(task->getChildNodes().size() == 2);
	mapStates_t::iterator itLeftChildCP = checkpoints.find(task->getChildNodes().front()->getId());
	assert(itLeftChildCP != checkpoints.end());
	mapStates_t::iterator itRightChildCP = checkpoints.find(task->getChildNodes().back()->getId());
	assert(itRightChildCP != checkpoints.end());

	// Draw new states for the chilrendre
	std::pair<size_t, size_t> childStates;

	if(ptrTensorCont->getOmega()->getDimensions()[0] == 0) { // no tensor: nothing to do
		childStates.first = parentState;
		childStates.second = parentState;
	} else {
		// Get the matrix corresponding to w_i
		const Tensor::eigenSparseMatrix_t &w_i = ptrTensorCont->getOmega()->getSparseTensor(event->getTime())[parentState];
		// Compute the product of this matrix and p_l(t) x p_r(t)
		Eigen::MatrixXd transitionMatrix = w_i.toDense();
		transitionMatrix = (transitionMatrix.array().colwise() * itLeftChildCP->second.second.getStateProb().array()).matrix();
		transitionMatrix = (transitionMatrix.array().rowwise() * itRightChildCP->second.second.getStateProb().array().transpose()).matrix();

		// Draw a state based on the derived probability transition matrix
		double randProb = transitionMatrix.sum() * Utils::RNG::Manager::getInstance()->getDefaultRNG().genUniformDbl();
		double sum = 0.;
		bool found = false;
		for(size_t iL=0; iL<(size_t)transitionMatrix.rows(); ++iL) {
			for(size_t iR=0; iR<(size_t)transitionMatrix.cols(); ++iR) {
				sum += transitionMatrix(iL, iR);
				if(sum >= randProb) {
					childStates.first = iL;
					childStates.second = iR;
					found = true;
					break;
				}
			}
			if(found) break;
		}
		assert(found && "Asynch speciation failed.");
	}

	// Set the children branches starting state
	// left children
	mapHistories_t::iterator itLeftChildHistory = histories.find(task->getChildNodes().front()->getEdgeIdMapping());
	assert(itLeftChildHistory == histories.end());
	histories[task->getChildNodes().front()->getEdgeIdMapping()].assign(1, std::make_pair(task->getTime(), childStates.first));
	// right children
	mapHistories_t::iterator itRightChildHistory = histories.find(task->getChildNodes().back()->getEdgeIdMapping());
	assert(itRightChildHistory == histories.end());
	histories[task->getChildNodes().back()->getEdgeIdMapping()].assign(1, std::make_pair(task->getTime(), childStates.second));
}

void SequentialStochasticMappingCPU::doSimulateEvent(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {

	Likelihood::Scheduler::Event* event = task->getEvent();

	if(event->checkEvent(Likelihood::Scheduler::PRESENT_TIME_EVENT)) {
		// Nothing to do here
	} else if(event->checkEvent(Likelihood::Scheduler::NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isSpeciationNode()) {
			// Recover parent edge last state
			assert(task->getParentNodes().size() == 1);
			mapHistories_t::iterator itParentHistory = histories.find(task->getParentNodes().front()->getEdgeIdMapping());
			assert(itParentHistory != histories.end());
			assert(!itParentHistory->second.empty());
			size_t parentState = itParentHistory->second.back().second;

			doSimulateAsyncSpeciationEvent(parentState, task, histories);
		} else if(eventNode->isSampledAncestor()) {
			// init the "downstream" segment with the upstream endState
			assert(task->getChildNodes().size() == 1 && task->getParentNodes().size() == 1);
			// Get parent (upstream) state:
			mapHistories_t::iterator itParentHistory = histories.find(task->getParentNodes().front()->getEdgeIdMapping());
			assert(itParentHistory != histories.end());
			assert(!itParentHistory->second.empty());
			size_t parentState = itParentHistory->second.back().second;

			// Set children (downstream) state:
			mapHistories_t::iterator itChildrenHistory = histories.find(task->getChildNodes().front()->getEdgeIdMapping());
			assert(itChildrenHistory == histories.end());
			histories[task->getChildNodes().front()->getEdgeIdMapping()].assign(1, std::make_pair(task->getTime(), parentState));
 		} else if(eventNode->isExtinct()) {
			// Nothing to do here
		}
	} else if(event->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT)) {

		// Get root prob
		Eigen::VectorXd observedRootFreq = probState.getStateProb().cwiseProduct(priorStateProbability);
		observedRootFreq /= observedRootFreq.sum();

		// Draw state
		size_t rootState = Utils::RNG::Manager::getInstance()->getDefaultRNG().getMultinomial(observedRootFreq);

		// Do speciation if needed
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isRootNode() && eventNode->isSpeciationNode()) {
			histories[eventNode->getId()].assign(2, std::make_pair(task->getTime(), rootState));
			doSimulateAsyncSpeciationEvent(rootState, task, histories);
		} else {
			// Set in history
			assert(task->getChildNodes().size() == 1);
			histories[task->getChildNodes().front()->getEdgeIdMapping()].assign(1, std::make_pair(task->getTime(), rootState));
		}
	} else if(event->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		// If its a speciation event
		if(eventNode->isSpeciationNode()) {
			// Recover parent edge last state
			assert(task->getParentNodes().size() == 1);
			mapHistories_t::iterator itParentHistory = histories.find(task->getParentNodes().front()->getEdgeIdMapping());
			assert(itParentHistory != histories.end());
			assert(!itParentHistory->second.empty());
			size_t parentState = itParentHistory->second.back().second;

			doSimulateAsyncSpeciationEvent(parentState, task, histories);
		} // else do nothing since we are still on the same branch
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT)) {
		// 1) lineage goes extinct: nothing to do
		// 2) lineage survives: nothing to do
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		// 1) lineage is sampled: new branch
		if(eventNode->isSampledAncestor()) {
			// init the "downstream" segment with the upstream endState
			assert(task->getChildNodes().size() == 1 && task->getParentNodes().size() == 1);
			// Get parent (upstream) state:
			mapHistories_t::iterator itParentHistory = histories.find(task->getParentNodes().front()->getEdgeIdMapping());
			assert(itParentHistory != histories.end());
			assert(!itParentHistory->second.empty());
			size_t parentState = itParentHistory->second.back().second;

			// Set children (downstream) state:
			mapHistories_t::iterator itChildrenHistory = histories.find(task->getChildNodes().front()->getEdgeIdMapping());
			assert(itChildrenHistory == histories.end());
			histories[task->getChildNodes().front()->getEdgeIdMapping()].assign(1, std::make_pair(task->getTime(), parentState));
		} // else: 2) lineage isn't sampled: nothing to do
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		// 1) lineage goes is sampled and goes extinct: nothing to do
		// 2) lineage isn't sampled: nothing to do
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RATE_SHIFT)) {
		// Nothing to do here
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RESCALING_EVENT)) {
		// Nothing to do here
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		// Nothing to do here
	} else {
		assert(false && "Event is not implemented.");
	}

	if(Utils::Output::outputManager().check(Utils::Output::HIGH_VERB) && !event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		std::cout << "Event :" << event->toString() << std::endl;
		std::cout << "State : " << probState.toString();
		std::cout << "----------------------------------------" << std::endl;
	}

}

void SequentialStochasticMappingCPU::doSimulateBranch(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {

	if(stochasticMappingAlgo == REJECTION_SAMPLING_ALGO) {
		bool succeeded = doRejectionSamplingAlgo(task, histories);
		// Rejection sampling failed
		if(!succeeded) {
			doDenseDopriAlgo(task, histories);
		}
	} else if(stochasticMappingAlgo == DENSE_EULER_ALGO) {
		doDenseEulerAlgo(task, histories);
	} else if(stochasticMappingAlgo == DENSE_DOPRI_ALGO) {
		doDenseDopriAlgo(task, histories);
	} else {
		assert(false && "Algorithm unknown.");
	}
}


bool SequentialStochasticMappingCPU::doDrawAncestralState(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {

	size_t edgeId = task->getEdgeIdMapping();
	// Get the history
	mapHistories_t::iterator itHistory = histories.find(edgeId);
	assert(itHistory != histories.end());

	// Get start and end state backward probabilities vectors
	mapStates_t::iterator itCheckpoint = checkpoints.find(task->getId());
	assert(itCheckpoint != checkpoints.end());

	// 1) Recover previous state from itHistory
	assert(!itHistory->second.empty());
	size_t curState = itHistory->second.back().second;		// State
	double curTime = itHistory->second.back().first;

	// 2) Define segment end points
	double startTime = task->getStartTime();
	double endTime = task->getEndTime();
	/*if(fabs(curTime-endTime) > (std::numeric_limits<double>::epsilon() * fabs(curTime+endTime) * 5.0)) {
		std::cout << startTime << std::endl;
		std::cout << endTime << std::endl;
		double endTime2 = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
		std::cout << endTime2 << std::endl;
		std::cout << curTime << std::endl;
		std::cout << fabs(curTime-endTime) << std::endl;
		std::cout << fabs(endTime2-endTime) << std::endl;
		std::cout << std::numeric_limits<double>::epsilon() * fabs(curTime+endTime) * 10.0 << std::endl;
		std::cout << task->toString() << std::endl;
		std::cout << "-----------------------------------" << std::endl;
	}*/
 	assert(fabs(curTime-endTime) <= 1.e-10);

	if(endTime != ptrScheduler->getEvents().back()->getTime()) {
		endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	}
	curTime = endTime; // This is necessary to match what is done in the backward phase

	// 3) Solve forward proability along the segment : f(t)
	stateType_t fwdProbState = itCheckpoint->second.first; 			// Init to something that make sense
	fwdProbState.getStateProb().setZero(); 							// Reset to zero
	fwdProbState.getStateProb()(curState) = 1.0;					// Set current state
	Specialized::AdaptiveSegmentFW fwApprox(startTime, endTime, fwdProbState, ptrTensorCont, ptrSApproxManager);

	// 4) Draw a state according to the "conditional" posterior distribution
	Eigen::VectorXd subtreeCondProb = itCheckpoint->second.first.getStateProb();
	Eigen::VectorXd posteriorProb = subtreeCondProb.cwiseProduct(fwApprox.getStateProb().getStateProb());
	posteriorProb /= posteriorProb.sum();
	size_t endState = Utils::RNG::Manager::getInstance()->getDefaultRNG().getMultinomial(posteriorProb);

	itHistory->second.push_back(std::make_pair(startTime, endState));

	return true;
}

bool SequentialStochasticMappingCPU::tryDrawingASample(double startTimeFW, double endTimeFW,
		size_t startStateFW, size_t endStateFW, mapHistories_t::iterator &itHistory) {

	// 1) Recovering U(t)
	Specialized::DenseUnobservedSharedPtr ptrDenseU = ptrSApproxManager->getDenseUnobserved(Specialized::FULL_PROCESS);

	// 2) Recover rates: rates can not change before the next DAG node
	Eigen::MatrixXd &lambda = ptrTensorCont->getEigenVecLambda(startTimeFW);
	Eigen::MatrixXd &mu = ptrTensorCont->getEigenVecMu(startTimeFW);
	Eigen::MatrixXd &phi = ptrTensorCont->getEigenVecPhi(startTimeFW);
	Eigen::MatrixXd &delta = ptrTensorCont->getEigenVecDelta(startTimeFW);
	Eigen::MatrixXd &eta = ptrTensorCont->getEigenMatrixEta(startTimeFW);

	// 3) Draw samples
	size_t cnt = 0;
	bool historyAccepted = false;
	while(!historyAccepted && cnt < REJECTION_SAMPLING_TRY_LIMIT) {
		 // reset to current time (oldest) and state
		double curTime = startTimeFW;
		size_t curState = startStateFW;

		// Draw a sample
		double acceptanceProb = 1.0;
		mapHistoriesVal_t tmpHistory;
		while(curTime > endTimeFW) {
			// Define rate: lambda_i + |eta_i,i|
			double rate = lambda(curState) - eta(curState, curState);

			// Draw waiting time
			double waitingTime = Utils::RNG::Manager::getInstance()->getDefaultRNG().genExponential(rate);
			waitingTime = std::min(curTime-endTimeFW, waitingTime);
			curTime -= waitingTime; // Decrease time

			// Condition on not observing ext./sampling events on the current sub-segment
			double sumUnobservedEventRates = mu(curState)+phi(curState)+delta(curState);
			if(sumUnobservedEventRates > 0.) {
				boost::math::poisson_distribution<double> poissonDistr(sumUnobservedEventRates*waitingTime);
				acceptanceProb *= boost::math::pdf(poissonDistr, 0);
			}

			// We reached the end of the segment
			if(curTime <= endTimeFW) {
				break; // or continue - same stuff
			}

			// Drawing an event
			double drawUnifDbl = Utils::RNG::Manager::getInstance()->getDefaultRNG().genUniformDbl();
			if(drawUnifDbl < lambda(curState)/rate) {
				// Speciation: we condition on having an unobserved sibling subtree out of the two possible branching outcomes
				Eigen::VectorXd u;
				ptrDenseU->defineProbabilities(curTime, u);
				double curU = u(curState);
				acceptanceProb *= 2.*curU;
			} else {
				// Transition: we draw a transition and new state according to x~eta_i,.
				Eigen::VectorXd transitionProbs = eta.row(curState);
				transitionProbs(curState) = 0.;
				transitionProbs /= transitionProbs.sum();
				curState = Utils::RNG::Manager::getInstance()->getDefaultRNG().getMultinomial(transitionProbs);
				std::pair<double, size_t> transition = std::make_pair(curTime, curState);
				tmpHistory.push_back(transition);
			}
		}

		// Rejection-sampling step
		double drawUnifDbl = Utils::RNG::Manager::getInstance()->getDefaultRNG().genUniformDbl();
		// If we end up in the right state and we don't reject the sample, then:
		if(curState == endStateFW && drawUnifDbl < acceptanceProb) {
			historyAccepted = true;
			// Add final state
			std::pair<double, size_t> finalTransition = std::make_pair(endTimeFW, curState);
			tmpHistory.push_back(finalTransition);
			// we register the history
			itHistory->second.insert(itHistory->second.end(), tmpHistory.begin(), tmpHistory.end());
			return true;
		}
		cnt++;
	}

	return false;
}

bool SequentialStochasticMappingCPU::doRejectionSamplingAlgo(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {

	// Time goes from oldest to youngest in this process => dt < 0

	// 1) Init
	// Recovering current history and checkpoints
	size_t edgeId = task->getEdgeIdMapping();
	mapHistories_t::iterator itHistory = histories.find(edgeId);
	assert(itHistory != histories.end());

	mapStates_t::iterator itCheckpoint = checkpoints.find(task->getId());
	assert(itCheckpoint != checkpoints.end());

	// Recover current state and time
	assert(!itHistory->second.empty());
	size_t curState = itHistory->second.back().second;	// Current state
	double curTime = itHistory->second.back().first;	// Current time (oldest)

	// 2 Define segment end points
	double startTimeBW = task->getStartTime();	// Target time (youngest)
	double endTimeBW   = task->getEndTime();	// Current time (oldest)
//	assert(fabs(curTime - endTimeBW) <= (std::numeric_limits<double>::epsilon() * fabs(curTime+endTimeBW) * 5.0));
	if(fabs(curTime - endTimeBW) <= (std::numeric_limits<double>::epsilon() * fabs(curTime+endTimeBW) * 5.0)) {
//		if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB) ) {
//			std::cout << "Warning: rejection sampling failed because the branch was too short." << std::endl;
//		}
		return false;
	}
	if(endTimeBW != ptrScheduler->getEvents().back()->getTime()) {
		endTimeBW = std::nextafter(endTimeBW,-std::numeric_limits<double>::infinity());
	}
	curTime = endTimeBW; // Required to match backward phase

	// 3) Solve forward probability along the segment
	stateType_t fwdProbState = itCheckpoint->second.first; 			// Init to something that make sense
	fwdProbState.getStateProb().setZero(); 							// Reset to zero
	fwdProbState.getStateProb()(curState) = 1.0;					// Set current state
	Specialized::AdaptiveSegmentFW fwApprox(startTimeBW, endTimeBW, fwdProbState, ptrTensorCont, ptrSApproxManager);

	// 4) Draw a state on the conditioned fw probability
	Eigen::VectorXd subtreeCondProb = itCheckpoint->second.first.getStateProb();
	Eigen::VectorXd posteriorProb = subtreeCondProb.cwiseProduct(fwApprox.getStateProb().getStateProb());
	posteriorProb /= posteriorProb.sum();
	size_t endState = Utils::RNG::Manager::getInstance()->getDefaultRNG().getMultinomial(posteriorProb);

	// 6) Check the expected rejection rate - if smaller that the "TRY_LIMIT" declare failure
	double expectedRejectionProb = fwApprox.getStateProb().getStateProb().dot(itCheckpoint->second.first.getStateProb() / itCheckpoint->second.first.getStateProb().sum());
	if(expectedRejectionProb < 1./(double)REJECTION_SAMPLING_TRY_LIMIT) {
		if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB) ) {
			std::cout << "Warning: rejection sampling failed because rejection probability = " << expectedRejectionProb << std::endl;
		}
		return false;
	}

	// Draw a sample endTimeBW -> startTimeBW
	bool succeeded = tryDrawingASample(endTimeBW, startTimeBW, curState, endState, itHistory);
	if(succeeded) {
		return true;
	} else {
		if(Utils::Output::outputManager().check(Utils::Output::MEDIUM_VERB)) {
			std::cout << "Warning: rejection sampling failed. No samples were found after " << REJECTION_SAMPLING_TRY_LIMIT << " tries." << std::endl;
		}
		// We rejected all samples => failure
		return false;
	}
}

bool SequentialStochasticMappingCPU::doDenseEulerAlgo(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {

	// Time goes from oldest to youngest in this process => dt < 0

	// 1) init
	size_t edgeId = task->getEdgeIdMapping();
	// Get the history
	mapHistories_t::iterator itHistory = histories.find(edgeId);
	assert(itHistory != histories.end());
	mapHistoriesVal_t tmpHistory;

	// Get start and end state backward probabilities vectors
	mapStates_t::iterator itCheckpoint = checkpoints.find(task->getId());
	assert(itCheckpoint != checkpoints.end());

	// 2) Recover segment times
	double startTimeBW = task->getStartTime();		// Youngest time
	double endTimeBW = task->getEndTime();			// Oldest time
	if(endTimeBW != ptrScheduler->getEvents().back()->getTime()) {
		endTimeBW = std::nextafter(endTimeBW,-std::numeric_limits<double>::infinity());
	}
    assert(fabs(itHistory->second.back().first - endTimeBW) <= (std::numeric_limits<double>::epsilon() * fabs(itHistory->second.back().first+endTimeBW) * 5.0));

	// 3) Recompute p(t) with a dense integrator for the current segment
	Specialized::DenseSegmentBW approxBwP(startTimeBW, endTimeBW, itCheckpoint->second.first, ptrTensorCont, ptrSApproxManager);
	// Check that the endpoint agrees with the checkpoint and the integator tolerance
	double absErr = (approxBwP.getStateProb().getStateProb() - itCheckpoint->second.second.getStateProb()).cwiseAbs().maxCoeff();
	if(absErr > 5.e-7) {
		std::cerr << "Absolute error = " << std::scientific << absErr << std::endl;
	}
	assert(absErr < 5.e-7);

	// 4) Forward simulate using an euler stepper
	euler_stepper_t eulerInt;

	// Init time, state and probabilities
	double curTime = endTimeBW;
	size_t curState = itHistory->second.back().second;
	stateType_t fwProbState = itCheckpoint->second.second;

	const double dt = FIXED_EULER_DT;
	while(curTime > startTimeBW) {

		// Reset the current state probability vector
		fwProbState.getStateProb().setZero();
		fwProbState.getStateProb()(curState) = 1.0;

		// Ensure we don't decrement over the end time (i.e., backward start time)
		double stepSize = std::min(dt, curTime-startTimeBW);

		// Do one Euler step
		eulerInt.do_step(intKernelFW, fwProbState, curTime, -stepSize);

		// Compute the conditional probability
		// compute p(t-dt)
		Eigen::VectorXd obsProb;
		approxBwP.defineProbabilities(curTime-stepSize, obsProb);
		// compute p(t-dt)*f(t-dt)
		Eigen::VectorXd transitionProbs = fwProbState.getStateProb().cwiseProduct(obsProb);
		transitionProbs /= transitionProbs.sum();

		// Draw new state
		size_t tmpNewState = Utils::RNG::Manager::getInstance()->getDefaultRNG().getMultinomial(transitionProbs);
		if(tmpNewState != curState) {
			curState = tmpNewState;
			// Save the transition
			std::pair<double, size_t> transition = std::make_pair(curTime, curState);
			tmpHistory.push_back(transition);
		}
		curTime -= dt;
	}

	// Save the end state and time
	std::pair<double, size_t> finalState = std::make_pair(startTimeBW, curState);
	tmpHistory.push_back(finalState);
	// Add the new events
	itHistory->second.insert(itHistory->second.end(), tmpHistory.begin(), tmpHistory.end());

	return true;
}

bool SequentialStochasticMappingCPU::doDenseDopriAlgo(Likelihood::Scheduler::DAG::NodeDAG* task, mapHistories_t &histories) {

	// Time goes from oldest to youngest in this process => dt < 0
	// 1) init
	size_t edgeId = task->getEdgeIdMapping();
	// Get the history
	mapHistories_t::iterator itHistory = histories.find(edgeId);
	assert(itHistory != histories.end());
	mapHistoriesVal_t tmpHistory;

	// Get start and end state backward probabilities vectors
	mapStates_t::iterator itCheckpoint = checkpoints.find(task->getId());
	assert(itCheckpoint != checkpoints.end());

	// 2) Recover segment times
	double startTimeBW = task->getStartTime();		// Youngest time
	double endTimeBW = task->getEndTime();			// Oldest time
	if(endTimeBW != ptrScheduler->getEvents().back()->getTime()) {
		endTimeBW = std::nextafter(endTimeBW,-std::numeric_limits<double>::infinity());
	}
    assert(fabs(itHistory->second.back().first - endTimeBW) <= (std::numeric_limits<double>::epsilon() * fabs(itHistory->second.back().first+endTimeBW) * 5.0));

	// 3) Recompute p(t) with a dense integrator for the current segment
	Specialized::DenseSegmentBW approxBwP(startTimeBW, endTimeBW, itCheckpoint->second.first, ptrTensorCont, ptrSApproxManager);
	// Check that the endpoint agrees with the checkpoint and the integator tolerance
	double absErr = (approxBwP.getStateProb().getStateProb() - itCheckpoint->second.second.getStateProb()).cwiseAbs().maxCoeff();
	if(absErr > 5.e-7) {
		std::cerr << "Absolute error = " << std::scientific << absErr << std::endl;
	}
	assert(absErr < 5.e-7);

	// Init time, state and probabilities
	double curTime = endTimeBW;
	size_t curState = itHistory->second.back().second;
	stateType_t fwProbState = itCheckpoint->second.second;	// Safe init for memory reasons
	stateType_t interpolatedFwProbState = fwProbState;		// Safe init for memory reasons

	while(curTime > startTimeBW) {

		// Get average dt for previous steps
		double avgAdaptiveDT = boost::accumulators::mean(accAdaptiveDT);
		avgAdaptiveDT = std::min(avgAdaptiveDT, curTime-startTimeBW);
		// Single step DOPRI integrator for f(t)
		Specialized::SingleStepSegmentFW ssSegFW(avgAdaptiveDT, ptrTensorCont, ptrSApproxManager);

		// reset f(t) to the current state
		fwProbState.getStateProb().setZero();
		fwProbState.getStateProb()(curState) = 1.0;

		// Do one adaptive step with a "dense" DOPRI integartor for f(t)
		double nextTime = ssSegFW.doNextStep(fwProbState, curTime, startTimeBW);
		// Adapt dt given current step
		accAdaptiveDT(nextTime-curTime);

		// Binary search
		const double TARGET_PROB = 1.-PROB_THRESHOLD;
		double stateProbBefore = 1.0;
		double timeBefore = curTime;
		double stateProbAfter = fwProbState.getStateProb()(curState);
		double timeAfter = nextTime;

		// Reduce one boundary of the time domain until we reach a value within then PROB_TOELRANCE
		while(fabs(stateProbBefore-TARGET_PROB) > PROB_TOLERANCE && fabs(stateProbAfter-TARGET_PROB) > PROB_TOLERANCE && stateProbAfter < TARGET_PROB) {
			// Linear approximation for the split time
			double splitTime = timeBefore - (timeBefore-timeAfter)*(stateProbBefore-TARGET_PROB)/(stateProbBefore-stateProbAfter);
			// Interpolate f(t) at splitTime
			ssSegFW.defineProbabilities(splitTime, interpolatedFwProbState);
			// Move the appropriate boundary to splitTime
			if(interpolatedFwProbState.getStateProb()(curState) < TARGET_PROB) {
				timeAfter = splitTime;
				stateProbAfter = interpolatedFwProbState.getStateProb()(curState);
			} else {
				timeBefore = splitTime;
				stateProbBefore = interpolatedFwProbState.getStateProb()(curState);
			}
		}

		double foundTime;
		// By default we always take the time further away from "curTime" and closer to "nextTime"
		if(stateProbAfter > TARGET_PROB || timeAfter == nextTime) {
			// We reached the end time or the rightmost boundary is higher than the target thershold probability
			// We use the adaptive integrator end point
			foundTime = timeAfter;
			interpolatedFwProbState = fwProbState; // copy th
		} else if(fabs(stateProbAfter-TARGET_PROB) < PROB_TOLERANCE) {
			// We use the rightmost endpoint of the binary search
			foundTime = timeAfter;
		} else {
			// We use the leftmost endpoint of the binary search
			foundTime = timeBefore;
		}

		// Define the posterior probability
		Eigen::VectorXd posteriorProb;
		// Get p(t)
		approxBwP.defineProbabilities(foundTime, posteriorProb);
		// Compute p(t)xf(t)
		Eigen::VectorXd transitionProbs = interpolatedFwProbState.getStateProb().cwiseProduct(posteriorProb);
		transitionProbs /= transitionProbs.sum();

		// Draw new state
		size_t tmpNewState = Utils::RNG::Manager::getInstance()->getDefaultRNG().getMultinomial(transitionProbs);
		if(tmpNewState != curState) {
			 // Transition, we reset the current state to the transition
			curState = tmpNewState;
			// Save the state
			std::pair<double, size_t> transition = std::make_pair(curTime, curState);
			tmpHistory.push_back(transition);
		}

		// Check if we reached the endpoint
		if(foundTime-startTimeBW <= 2.*std::numeric_limits<double>::epsilon()) {
			// Reset to endpoint time: required for numerical accuracy
			curTime = startTimeBW;
		} else {
			// Set current time to the new time
			curTime = foundTime;
		}
	}

	// Save last state
	std::pair<double, size_t> finalState = std::make_pair(startTimeBW, curState);
	tmpHistory.push_back(finalState);

	// Augment the stochastic map with the new events
	itHistory->second.insert(itHistory->second.end(), tmpHistory.begin(), tmpHistory.end());

	return true;
}

}
} /* namespace Approximator */
} /* namespace Likelihood */
