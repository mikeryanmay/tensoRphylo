/*
 * ParallelSchedulerDAG.cpp
 *
 *  Created on: Apr 14, 2020
 *      Author: xaviermeyer
 */

#include "ParallelSchedulerDAG.h"

#include <vector>

#include "../../../Data/Structure/Node.h"
#include "../Event.h"
#include "Utils/Parallel/Manager.h"
#if defined(_OPENMP)

#include <set>
#include <iostream>

#include "Likelihood/StateTypes/Vector/EigenState.h"
#include "NodeDAG.h"
#include "DAG.h"

namespace Likelihood {
namespace Scheduler {
namespace DAG {

ParallelSchedulerDAG::ParallelSchedulerDAG(size_t nThread, DAGSharedPtr aPtrDAG) : N_THREADS(nThread), ptrDAG(aPtrDAG) {
	assert(N_THREADS >= 1);
	init();
}

ParallelSchedulerDAG::~ParallelSchedulerDAG() {
	for(size_t iL=0; iL<locks.size(); ++iL) {
		omp_destroy_lock(&locks[iL]);
	}
}

bool ParallelSchedulerDAG::isDone() {
	return rootSignaled;
}

void ParallelSchedulerDAG::reset() {

	assert(ptrDAG->getRoot()->getProbStates().size() == 0);
	init();
}

NodeDAG* ParallelSchedulerDAG::getNextTask(size_t iThread) {
	assert(iThread < N_THREADS);

	NodeDAG* task = getNextTaskFromLocalThreadPool(iThread);
	//std::cout << "[ " << iThread << "] Trying to get a task from local thread pool : " << (task == NULL ? "failure" : "success") << std::endl;
	if(task == NULL) {
		task = getNextTaskFromGlobalThreadPool(iThread);
		//std::cout <<  "[ " << iThread << "] Trying to get a task from global thread pool : " << (task == NULL ? "failure" : "success") << std::endl;
	}

	return task;
}

NodeDAG* ParallelSchedulerDAG::getNextTaskFromLocalThreadPool(size_t iThread) {
	assert(N_THREADS >= 1 && iThread < N_THREADS);

	NodeDAG* task = NULL;
	omp_set_lock(&locks[iThread]);

	if(!pendingTasks[iThread].empty()) {
		task = pendingTasks[iThread].front();
		pendingTasks[iThread].pop_front();
	}

	omp_unset_lock(&locks[iThread]);

	return task;
}

NodeDAG* ParallelSchedulerDAG::getNextTaskFromGlobalThreadPool(size_t iThread) {
	assert(N_THREADS >= 1 && iThread < N_THREADS);

	assert(pendingTasks[iThread].empty());

	for(size_t iP=0; iP<N_THREADS; ++iP) {
		if(iP == iThread) continue;

		NodeDAG* task = NULL;
		omp_set_lock(&locks[iP]);
		if(!pendingTasks[iP].empty()) {
			task = pendingTasks[iP].front();
			pendingTasks[iP].pop_front();
		}
		omp_unset_lock(&locks[iP]);

		if(task != NULL) return task;
	}

	return NULL;
}

bool ParallelSchedulerDAG::tryNonBlockingSignalAndPop(size_t iThread, NodeDAG* &inOutNode) {
	// If it's the root, we are done
	if(inOutNode == ptrDAG->getRoot()) {
		rootSignaled = true;
		ptrDAG->getRoot()->resetChildReadyState();
		inOutNode = NULL;
		return true;
	}

	// Else we check who can be added to the pool of pending task
	// The following can only be one of both for now but if dense integrator would be added to the dag: the following wouldn't work
	NodeDAG *aNode = inOutNode;

	NodeDAG *parentNode = aNode->getParentNodes().front();

	bool isParentReady = false;
	#pragma omp critical (altering_nodes)
	{
		assert((aNode->isEventTask() || aNode->isIntegrationTask()) && aNode->getParentNodes().size() == 1 && aNode->getProbStates().size() == 1);
		aNode->getProbStates().front()->setEdgeMapping(parentNode->getEdgeIdMapping());
		parentNode->attachState(aNode->getProbStates().front());

		aNode->detachStates();

		aNode->signalReadyToParents();
		aNode->resetChildReadyState();

		isParentReady = parentNode->isReady();
	}

	if(isParentReady) {
		inOutNode = parentNode;
		return true;
	} else {
		inOutNode = NULL;
		return false;
	}

}


size_t ParallelSchedulerDAG::getNThread() const {
	return N_THREADS;
}

void ParallelSchedulerDAG::signalDone(size_t iThread, NodeDAG* aNode) {
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

	NodeDAG *parentNode = aNode->getParentNodes().front();

	bool isParentReady = false;
	#pragma omp critical (altering_nodes)
	{
		assert((aNode->isEventTask() || aNode->isIntegrationTask()) && aNode->getParentNodes().size() == 1 && aNode->getProbStates().size() == 1);

		aNode->getProbStates().front()->setEdgeMapping(parentNode->getEdgeIdMapping());
		parentNode->attachState(aNode->getProbStates().front());

		aNode->detachStates();

		aNode->signalReadyToParents();
		aNode->resetChildReadyState();

		isParentReady = parentNode->isReady();
	}

	if(isParentReady) {
		omp_set_lock(&locks[iThread]);
		pendingTasks[iThread].push_back(parentNode);
		omp_unset_lock(&locks[iThread]);
	}
}

void ParallelSchedulerDAG::init() {

	rootSignaled = false;
	pendingTasks.clear();
	pendingTasks.resize(N_THREADS);

	if(locks.empty()) {
		locks.resize(N_THREADS);
		for(size_t iL=0; iL<locks.size(); ++iL) {
			omp_init_lock(&locks[iL]);
		}
	}

	{
		// Find N_THREADS sub roots to split the initial work
		std::list<NodeDAG*> subtreesRoot(1, ptrDAG->getRoot());
		while(subtreesRoot.size() < N_THREADS && subtreesRoot.size() < ptrDAG->getTips().size()) {
			// get and clean front
			NodeDAG* topNode = subtreesRoot.front();
			subtreesRoot.pop_front();

			if(topNode->getChildNodes().empty()) {
				subtreesRoot.push_back(topNode);
			} else {
				subtreesRoot.insert(subtreesRoot.begin(), topNode->getChildNodes().begin(), topNode->getChildNodes().end());
			}
		}

		// Subroot are found: get the starting nodes
		size_t iThread = 0;
		std::set<size_t> visitedTips; // We need to kip track of that since it's a DAG
		for(std::list<NodeDAG*>::iterator itR = subtreesRoot.begin(); itR != subtreesRoot.end();  ++itR) {
			std::vector<NodeDAG*> subsetTips;
			(*itR)->identifyTips(subsetTips, visitedTips);

			pendingTasks[iThread].insert(pendingTasks[iThread].end(), subsetTips.begin(), subsetTips.end());
			iThread++;
		}
	}

	// Find N_THREADS sub roots to split the initial work
	/*std::list<std::pair< NodeDAG*, double > > subtreesRoot(1, std::make_pair(ptrDAG->getRoot(), ptrDAG->getRoot()->identifySubDAGSumBL()));
	//std::cout << "Total BL : " << ptrDAG->getRoot()->identifySubDAGSumBL() << std::endl;
	while(subtreesRoot.size() < 4*N_THREADS && subtreesRoot.size() < ptrDAG->getTips().size()) {
		// get and clean front
		NodeDAG* topNode = subtreesRoot.front().first;
		subtreesRoot.pop_front();

		if(topNode->getChildNodes().empty()) {
			subtreesRoot.push_back(std::make_pair(topNode, topNode->identifySubDAGSumBL()));
		} else {
			for(size_t iC=0; iC<topNode->getChildNodes().size(); ++iC) {
				subtreesRoot.push_back(std::make_pair(topNode->getChildNodes()[iC], topNode->getChildNodes()[iC]->identifySubDAGSumBL()));
			}
		}
	}


	std::vector<double> sumBLPerThread(N_THREADS, 0.);
	std::set<size_t> visitedTips; // We need to kip track of that since it's a DAG
	for(std::list< std::pair< NodeDAG*, double > >::iterator itR = subtreesRoot.begin(); itR != subtreesRoot.end();  ++itR) {
		size_t iMinT = 0;
		size_t iMaxT = 0;
		for(size_t iT=1; iT<N_THREADS; ++iT) {
			if(sumBLPerThread[iT] < sumBLPerThread[iMinT]) {
				iMinT = iT;
			}
			if(sumBLPerThread[iT] > sumBLPerThread[iMaxT]) {
				iMaxT = iT;
			}
		}

		std::vector<NodeDAG*> subsetTips;
		(*itR).first->identifyTips(subsetTips, visitedTips);

		pendingTasks[iMinT].insert(pendingTasks[iMinT].end(), subsetTips.begin(), subsetTips.end());
		sumBLPerThread[iMinT] += (*itR).second;
	}*/

	/*for(size_t i=0; i<N_THREADS; ++i) {
		std::cout << "Work thread " << i << " :: " << sumBLPerThread[i] << std::endl;
	}
	std::cout << "--------------------------" << std::endl;*/


}

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */

#endif
