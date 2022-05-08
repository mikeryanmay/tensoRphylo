/*
 * SequentialOptimizedCPU.cpp
 *
 *  Created on: Nov 19, 2019
 *      Author: xaviermeyer
 */

#include "SequentialOptimizedCPU.h"

#include <boost/core/ref.hpp>
#include <boost/numeric/odeint/util/unit_helper.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <cmath>

#include "Likelihood/CustomIntegrators/AdaptiveIntegrators.hpp"
#include "Utils/Output/OutputManager.h"
#include "Tensor/IncTensor.h"
#include "SynchronousEvents/IncSynchronousEvents.h"
#include "Data/Reader/IncPhyloReader.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Likelihood/Scheduler/IncScheduler.h"

namespace Likelihood {
namespace Approximator {

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::SequentialOptimizedCPU(
								Likelihood::Integrator::integrationScheme_t aIntScheme,
								Conditions::conditionalProbability_t aConditionType,
								Phylogeny::Data::ContainerSharedPtr aPtrData,
								Scheduler::SchedulerSharedPtr aPtrScheduler,
								SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
								Tensor::ContainerSharedPtr aPtrTensorCont) :
								SynchronousApproximator(aIntScheme, aConditionType, aPtrData, aPtrScheduler, aPtrSyncEventsCont, aPtrTensorCont),
										N_MAX_STATE_VECTOR(ptrScheduler->getPtrTree()->getTerminalNodes().size()+Likelihood::Conditions::getNVector(conditionalProbType)),
										kernels(ptrData, ptrSyncEventsCont, ptrTensorCont, condCompatibilityMode),
										intKernel(N_MAX_STATE_VECTOR, ptrData, ptrSyncEventsCont, ptrTensorCont) {

	assert((!withCladoEvents  || ptrTensorCont->getOmega()->isSparse()) && "This kernel only works with sparse cladogenetic matrices since commit '6296fe4'.");

	// Important: initialize the size of vector that the memory pull will allocate.
	Utils::MemoryPool::eigenCPU().setNCategories(ptrTensorCont->getNumberOfState());
	Utils::MemoryPool::eigenCPU().setMaxStatesVector(N_MAX_STATE_VECTOR);

	ptrIntegrator = Likelihood::Integrator::Factory::createIntegrator<stateType_t, intKernel_t, operations_t>(DEFAULT_ABS_TOLERANCE, DEFAULT_REL_TOLERANCE, deltaT, aIntScheme);

}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::~SequentialOptimizedCPU() {
	_REPORT_TO_FILE(Utils::Profiling::UniqueProfiler::getInstance())
	assert(ptrIntegrator);
	delete ptrIntegrator;
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::setDefaultDeltaT(double aDeltaT) {
	assert(ptrIntegrator != NULL);
	deltaT = aDeltaT;
	ptrIntegrator->setDeltaT(deltaT);
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
size_t SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::getTotalNumberOfIntegrationSteps() const {
	assert(ptrIntegrator != NULL);
	return ptrIntegrator->getNSteps();
}


template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::doPreProcessingSteps() {
	// Init with an empty state
	//Likelihood::Integrator::RungeKuttaDOPRI5<stateType_t, intKernel_t, operations_t> tmp(DEFAULT_ABS_TOLERANCE, DEFAULT_REL_TOLERANCE, 1.0);
	//integratorRKD5 = tmp;

	probState = Likelihood::StateType::Optimized::EigenState<conditionalProbType>();
	integrationTimes.clear();

	vecProbesState.clear();

	/*std::cout << "Reinit probState : " << std::endl;
	std::cout << probState.toString() << std::endl;
	std::cout << "----------------------------------------" << std::endl;*/
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::doIntegrationStep(size_t iEdgesLayer) {
	// We can get the layer of edges using the scheduler
	// However, it's not even needed as the kernels already deal with that
	// After a doEventStep, the state has already been resized to fit the next set of edges

	// iEdgesLayer starts at iEvent = iEdgesLayer and ends at iEvent+1=iEdgesLayer+1
	double startTime = ptrScheduler->getEvents()[iEdgesLayer]->getTime();
	double endTime = ptrScheduler->getEvents()[iEdgesLayer+1]->getTime();

	if(endTime == startTime) return;

	if(endTime != ptrScheduler->getEvents().back()->getTime()) {
		endTime = std::nextafter(endTime,-std::numeric_limits<double>::infinity());
	}
	//std::cout << "Staring event : " << ptrScheduler->getEvents()[iEdgesLayer]->toString() << std::endl;
	//std::cout << "Ending event : " << ptrScheduler->getEvents()[iEdgesLayer+1]->toString() << std::endl;

	assert(ptrIntegrator != NULL);
	ptrIntegrator->integrate(startTime, endTime, probState, intKernel);

	integrationTimes.insert(integrationTimes.end(), ptrIntegrator->getVecTimes().begin(), ptrIntegrator->getVecTimes().end());

	/*std::cout << "After integration : " << std::endl;
	std::cout << probState.toString() << std::endl;
	std::cout << "----------------------------------------" << std::endl;*/


}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::doEventStep(size_t iEvent) {

	/*std::cout << "Before event : " << iEvent << std::endl;
	std::cout << probState.toString() << std::endl;
	std::cout << "----------------------------------------" << std::endl;*/

	Likelihood::Scheduler::Event* event = ptrScheduler->getEvents()[iEvent];

	if(event->checkEvent(Likelihood::Scheduler::PRESENT_TIME_EVENT)) {			// When we get a sampled ancestor, we need to be careful when we are in the single branch hack case
		kernels.setInitialCondition(event->getNodes(), probState);
	} else if(event->checkEvent(Likelihood::Scheduler::NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isSpeciationNode()) {
			kernels.computeAsynchSpeciation(event->getTime(), eventNode, probState);
		} else if(eventNode->isSampledAncestor()) {
			// When we get a sampled ancestor, we need to be careful when we are in the single branch hack case
			kernels.computeAsynchSampling(event->getTime(), eventNode, probState);
			kernels.rescaleRequestedBranchesProbabilities(event->getNodes(), probState); // rescale here
		} else if(eventNode->isExtinct()) {
			kernels.setInitalExtinctNodeCondition(event->getTime(), eventNode, probState);
		}
	} else if(event->checkEvent(Likelihood::Scheduler::FINAL_NODE_EVENT)) {
		PS::Node* eventNode = event->getNodes()[0];
		if(eventNode->isRootNode() && eventNode->isSpeciationNode()) {
			kernels.computeAsynchSpeciation(event->getTime(), event->getNodes().front(), probState);
		}
		Eigen::VectorXd rf = this->getRootFrequency(event->getTime());
		logLikelihood = kernels.computeLogLikelihood(event->getTime(), eventNode, rf, probState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYCHRONOUS_SPECIATION_EVENT)) {
		kernels.computeMassSpeciation(event->getTime(), event->getNodes(), probState);
		probState.rescaleAll();
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_EXTINCTION_EVENT)) {
		kernels.computeMassExtinction(event->getTime(), event->getNodes(), probState);
		probState.rescaleAll();
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_SAMPLING_EVENT)) {
		kernels.computeMassSamplingEvent(event->getTime(), event->getNodes(), probState);
		probState.rescaleAll();
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		kernels.computeMassDestrSamplingEvent(event->getTime(), event->getNodes(), probState);
		probState.rescaleAll();
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RATE_SHIFT)) {
		// DO NOTHING BUT ENSURE THAT THE NUMERICAL INTEGRATOR DEAL ONLY WITH CONTIUNOUS FUNCTIONS
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_RESCALING_EVENT)) {
		kernels.rescaleRequestedBranchesProbabilities(event->getNodes(), probState);
	} else if(event->checkEvent(Likelihood::Scheduler::SYNCHRONOUS_MONITORING_PROB)) {
		this->doReportState(event->getTime());
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

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::doPostProcessingSteps() {
	// Nothing to be done?
}

template < Likelihood::Conditions::conditionalProbability_t conditionalProbType, Tensor::etaStructure_t etaStructure, bool withCladoEvents >
void SequentialOptimizedCPU<conditionalProbType, etaStructure, withCladoEvents>::doReportState(double t) {

	Likelihood::Monitor::ProbeState probe;
	probe.time = t;
	probe.u = probState.getUnobservedStateProb();

	for(size_t i=0; i<probState.size(); ++i) {
		probe.vecP.push_back(probState.getObservedStateProb().col(i));
		probe.vecIdEdge.push_back(probState.getVecProbToEdgeMapping()[i]);
		probe.scalingFactor.push_back(probState.getScalingFactorByVecPos(i));
	}
	std::pair<double, bool> resSR = intKernel.getStiffnessRatio(t, probState.getUnobservedStateProb());
	probe.stiffnessRatio = resSR.first;
	probe.hasNegativeImgEIGVal = resSR.second;

	vecProbesState.push_back(probe);
}


} /* namespace Approximator */
} /* namespace Likelihood */
