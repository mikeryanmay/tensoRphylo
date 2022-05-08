/*
 * BranchwiseCPU.cpp
 *
 *  Created on: April 17, 2020
 *      Author: xaviermeyer
 */

#include "SequentialBranchwiseCPU.h"

#include <algorithm>
#include <cmath>

#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "IncLikelihoodApproximator.h"
#include "Likelihood/Scheduler/IncScheduler.h"
#include "Utils/Output/OutputManager.h"
#include "Likelihood/Scheduler/DAG/IncDAG.h"
#include "Likelihood/CustomIntegrators/AdaptiveIntegrators.hpp"

namespace Likelihood {
namespace Approximator {

SequentialBranchwiseCPU::SequentialBranchwiseCPU(Likelihood::Integrator::integrationScheme_t aIntScheme,
	   	   	   	   	   	   	    Conditions::conditionalProbability_t aConditionType,
								Phylogeny::Data::ContainerSharedPtr aPtrData,
								Scheduler::SchedulerSharedPtr aPtrScheduler,
								SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
								Tensor::ContainerSharedPtr aPtrTensorCont) :
		AsynchronousApproximator(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
		kernels(aConditionType, ptrData, ptrSApproxManager, ptrSyncEventsCont, ptrTensorCont, condCompatibilityMode),
		intKernel(ptrSApproxManager, ptrTensorCont) {

	// Important: initialize the size of vector that the memory pull will allocate.
	Utils::MemoryPool::eigenCPU().setNCategories(ptrTensorCont->getNumberOfState());
	Utils::MemoryPool::eigenCPU().setMaxStatesVector(ptrScheduler->getPtrTree()->getTerminalNodes().size());

	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<stateType_t, intKernel_t, operations_t >(DEFAULT_ABS_TOLERANCE, DEFAULT_REL_TOLERANCE, deltaT, aIntScheme);

}

SequentialBranchwiseCPU::~SequentialBranchwiseCPU() {
	delete ptrIntegrator;
}

void SequentialBranchwiseCPU::setDefaultDeltaT(double aDeltaT) {
	assert(ptrIntegrator != NULL);
	deltaT = aDeltaT;
	ptrIntegrator->setDeltaT(deltaT);
}

size_t SequentialBranchwiseCPU::getTotalNumberOfIntegrationSteps() const {
	return ptrIntegrator->getNSteps();
}

void SequentialBranchwiseCPU::doPreProcessingSteps() {
	// Init with an empty state
	probState = Likelihood::StateType::Vector::EigenState();

	if(ptrScheduler->hasBeenUpdated()) {
		// reset DAG
		ptrDAG->rebuild();
		ptrScheduler->clearHasBeenUpdatedFlag();
	}

	ptrSchedDAG->reset();

	integrationTimes.clear();

	vecProbesState.clear();

}

void SequentialBranchwiseCPU::doIntegrationStep(Likelihood::Scheduler::DAG::NodeDAG* task) {
	// We can get the layer of edges using the scheduler
	// However, it's not even needed as the kernels already deal with that
	// After a doEventStep, the state has already been resized to fit the next set of edges

	// iEdgesLayer starts at iEvent = iEdgesLayer and ends at iEvent+1=iEdgesLayer+1
	assert(task->isIntegrationTask());
	double startTime = task->getStartTime();
	double endTime = task->getEndTime();

	if(endTime == startTime) return;

	// We make sure that endTime is always within -1*epsilon from the true value
	// to make sure that we can find in which interval we are
	endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	if(endTime != ptrScheduler->getEvents().back()->getTime()) {
		endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	}


	/*std::cout << "Before integration : " << std::endl;
	std::cout << probState.toString() << std::endl;
	std::cout << "----------------------------------------" << std::endl;*/
	assert(task->getProbStates().size() == 1);
	stateType_t &activePState = *task->getProbStates().front();

	assert(ptrIntegrator != NULL);
	ptrIntegrator->integrate(startTime, endTime, activePState, intKernel);

	integrationTimes.insert(integrationTimes.end(), ptrIntegrator->getVecTimes().begin(), ptrIntegrator->getVecTimes().end());

	/*std::cout << "After integration : " << std::endl;
	std::cout << probState.toString() << std::endl;
	std::cout << "----------------------------------------" << std::endl;*/
}

void SequentialBranchwiseCPU::doEventStep(Likelihood::Scheduler::DAG::NodeDAG* task) {

	/*std::cout << "Before event : " << std::endl;
	std::cout << probState.toString() << std::endl;
	if(conditionType == Likelihood::Conditions::STEM_TWO_EXT_SAMPLES) {
		std::cout << "U_hat : " << probState.getUnobservedNoSamplingStateProb().transpose() << std::endl;
		std::cout << "S_hat : " << probState.getSingletonNoSamplingStateProb().transpose() << std::endl;
	}
	std::cout << "----------------------------------------" << std::endl;*/

	Likelihood::Scheduler::Event* event = task->getEvent();

	if(event->checkEvent(Likelihood::Scheduler::PRESENT_TIME_EVENT)) {			// When we get a sampled ancestor, we need to be careful when we are in the single branch hack case
		kernels.setInitialCondition(task);
	} else if(event->checkEvent(Likelihood::Scheduler::NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isSpeciationNode()) {
			kernels.computeAsynchSpeciation(event->getTime(), task);
		} else if(eventNode->isSampledAncestor()) {
			// When we get a sampled ancestor, we need to be careful when we are in the single branch hack case
			kernels.computeAsynchSampling(event->getTime(), task);
		} else if(eventNode->isExtinct()) {
			kernels.setInitalExtinctNodeCondition(event->getTime(), task);
		}
	} else if(event->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isRootNode() && eventNode->isSpeciationNode()) {
			kernels.computeAsynchSpeciation(event->getTime(), task);
		}
		Eigen::VectorXd rf = this->getRootFrequency(event->getTime());
		logLikelihood = kernels.computeLogLikelihood(event->getTime(), rf, task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT)) {
		kernels.computeMassSpeciation(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT)) {
		kernels.computeMassExtinction(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT)) {
		kernels.computeMassSamplingEvent(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		kernels.computeMassDestrSamplingEvent(event->getTime(), task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RATE_SHIFT)) {
		// DO NOTHING BUT ENSURE THAT THE NUMBERICAL INTEGRATOR DEAL ONLY WITH CONTIUNOUS FUNCTIONS
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RESCALING_EVENT)) {
		kernels.rescaleProbabilities(task);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		kernels.rescaleProbabilities(task); // Forcing rescaling before probe
		this->doReportState(task);
	} else {
		assert(false && "Event is not implemented.");
	}

	if(Utils::Output::outputManager().check(Utils::Output::HIGH_VERB) && !event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		std::cout << "Event :" << event->toString() << std::endl;
		std::cout << "State : " << probState.toString();
		std::cout << "----------------------------------------" << std::endl;
	}

	/*std::cout << "After event : " << std::endl;
	std::cout << probState.toString() << std::endl;
	std::cout << "----------------------------------------" << std::endl;*/

}

void SequentialBranchwiseCPU::doPostProcessingSteps() {
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

void SequentialBranchwiseCPU::doReportState(Likelihood::Scheduler::DAG::NodeDAG* task) {

	if(task->getProbStates().empty()) return;

	Likelihood::Monitor::ProbeState probe;
	probe.time = task->getTime();
	ptrSApproxManager->getDenseUnobserved(Likelihood::Approximator::Specialized::FULL_PROCESS)->defineProbabilities(probe.time, probe.u);

	for(size_t iP=0; iP<task->getProbStates().size(); ++iP) {
		probe.vecP.push_back(task->getProbStates()[iP]->getStateProb());
		probe.vecIdEdge.push_back(task->getProbStates()[iP]->getEdgeMapping());
		probe.scalingFactor.push_back(task->getProbStates()[iP]->getScaling());
	}

	std::pair<double, bool> resSR = intKernel.getStiffnessRatio(probe.time, probe.u);
	probe.stiffnessRatio = resSR.first;
	probe.hasNegativeImgEIGVal = resSR.second;

	vecProbesState.push_back(probe);
}


} /* namespace Approximator */
} /* namespace Likelihood */
