/*
 * SchedulerDAG.cpp
 *
 *  Created on: Apr 14, 2020
 *      Author: xaviermeyer
 */

#include "SchedulerDAG.h"

#include <set>
#include <iostream>

#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "NodeDAG.h"
#include "DAG.h"

namespace Likelihood {
namespace Scheduler {
namespace DAG {

SchedulerDAG::SchedulerDAG(DAGSharedPtr aPtrDAG) : ptrDAG(aPtrDAG) {
	init();
}

SchedulerDAG::~SchedulerDAG() {

}

bool SchedulerDAG::isDone() {
	return rootSignaled;
}

void SchedulerDAG::reset() {

	assert(ptrDAG->getRoot()->getProbStates().size() == 0);
	init();
}

NodeDAG* SchedulerDAG::getNextTask() {
	assert(!rootSignaled);

	NodeDAG* task = pendingTasks.front();
	pendingTasks.pop_front();

	return task;
}

void SchedulerDAG::signalDone(NodeDAG* aNode) {
	// If it's the root, we are done
	if(aNode == ptrDAG->getRoot()) {
		rootSignaled = true;
		ptrDAG->getRoot()->resetChildReadyState();
		return;
	}

	// Else we check who can be added to the pool of pending task
	// The following can only be one of both for now but if dense integrator would be added to the dag: the following wouldn't work
	//std::cout << ptrDAG->getRoot()->toString() << std::endl;
	//std::cout << aNode->toString() << std::endl;
	assert((aNode->isEventTask() || aNode->isIntegrationTask()) && aNode->getParentNodes().size() == 1 &&  aNode->getProbStates().size() == 1);

	NodeDAG *parentNode = aNode->getParentNodes().front();

	aNode->getProbStates().front()->setEdgeMapping(parentNode->getEdgeIdMapping());
	parentNode->attachState(aNode->getProbStates().front());

	aNode->detachStates();

	aNode->signalReadyToParents();
	aNode->resetChildReadyState();

	if(parentNode->isReady()) {
		pendingTasks.push_back(parentNode);
	}
}

void SchedulerDAG::init() {

	rootSignaled = false;
	pendingTasks.clear();

	// Just copy dag ready nodes
	pendingTasks.insert(pendingTasks.end(), ptrDAG->getTips().begin(), ptrDAG->getTips().end());

}

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */
