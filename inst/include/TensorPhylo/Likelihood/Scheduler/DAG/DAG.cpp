/*
 * DAG.cpp
 *
 *  Created on: Apr 13, 2020
 *      Author: xaviermeyer
 */

#include "DAG.h"

#include <cmath>
#include <cstdio>

#include "../../Approximator/BaseApproximator.h"
#include "../BaseScheduler.h"
#include "../Event.h"
#include "NodeDAG.h"
#include "Data/Structure/IncTreeStructure.h"

using Likelihood::Approximator::BaseApproximator;


namespace Likelihood {
namespace Scheduler {
namespace DAG {

DAG::DAG(Scheduler::SchedulerSharedPtr aPtrScheduler) : ptrScheduler(aPtrScheduler) {
	assert(ptrScheduler);

	root = NULL;
	build();
	ptrScheduler->clearHasBeenUpdatedFlag();
}

DAG::~DAG() {

	clear();

}

void DAG::initTips(Event* currentEvent, mappingEdgeToDAG_t &edgeToDAG) {

	// first init nodes:
	assert(currentEvent->checkEvent(PRESENT_TIME_EVENT));
	// Create tips nodes (one per branch
	for(size_t iN=0; iN<currentEvent->getNodes().size(); ++iN) {
		Event* branchEvent = new Event(PRESENT_TIME_EVENT, currentEvent->getNodes()[iN]);
		NodeDAG *node = new NodeDAG(currentEvent->getTime(), branchEvent, currentEvent->getNodes()[iN]->getEdgeToParent()->getId());
		tips.push_back(node);
		edgeToDAG[currentEvent->getNodes()[iN]->getEdgeToParent()->getId()] = node;
	}

}


NodeDAG* DAG::addSegmentedIntegrationNodes(double startTime, double endTime, size_t edgeId, NodeDAG *prevNode) {

	double edgeLength = endTime-startTime;
	if(edgeLength > BaseScheduler::MAX_SEGMENT_SIZE_WITHOUT_RESCALING) {
		//std::cout << "Fragmenting branch of length : " << edgeLength << " ( " << startTime << " -> " << endTime << " ) " << std::endl;
		size_t nSubSegment = std::floor(edgeLength/BaseScheduler::MAX_SEGMENT_SIZE_WITHOUT_RESCALING);
		double subsegmentSize = edgeLength/nSubSegment;

		// Request a rescaling at time child_node + i*subsegment_length
		double sTime = startTime;
		for(size_t iS=0; iS<nSubSegment; ++iS) {
			double eTime = std::min(endTime, sTime+subsegmentSize);

			//std::cout << "Fragment " << iS << " : " << (eTime-sTime) << " ( " << sTime << " -> " << eTime << " ) " << std::endl;
			// Create integration
			NodeDAG *newNode = new NodeDAG(sTime, eTime, edgeId);
			prevNode->addParentNode(newNode);
			newNode->addChildNode(prevNode);

			sTime = eTime;
			prevNode = newNode;

			// Create scaling node if needed
			if(iS < nSubSegment-1) {
				//std::cout << "Adding a rescaling event at time " << eTime << std::endl;
				Event* rescaleEvent = new Event(SYNCHRONOUS_RESCALING_EVENT, eTime);
				NodeDAG *newNode = new NodeDAG(eTime, rescaleEvent, edgeId);
				prevNode->addParentNode(newNode);
				newNode->addChildNode(prevNode);
				prevNode = newNode;
			}
		}
		//getchar();
	} else {
		//std::cout << "No fragment : " << (endTime-startTime) << std::endl;
		NodeDAG *newNode = new NodeDAG(startTime, endTime, edgeId);
		prevNode->addParentNode(newNode);
		newNode->addChildNode(prevNode);

		prevNode = newNode;
	}

	return prevNode;
}

void DAG::asynchronousEvents(Event* currentEvent, mappingEdgeToDAG_t &edgeToDAG) {
	assert(currentEvent->getNodes().size() == 1 && currentEvent->getNodes().front() != NULL);
	PS::Node* eventNodeTree = currentEvent->getNodes().front();

	//std::cout << "Async event : " << currentEvent->toString() << std::endl;
	//std::cout << "Async event [N] : " << eventNodeTree->toString() << std::endl;

	// Creating the dag node and mapping with previous DAG node if existing: add integration step too
	Event* branchEvent = new Event(currentEvent->getEvent(), currentEvent->getNodes().front());

	long int branchId = -1;
	if(!eventNodeTree->isOriginNode()) {
		branchId = currentEvent->getNodes().front()->getEdgeToParent()->getId();
	}
	NodeDAG *eventNodeDAG = new NodeDAG(currentEvent->getTime(), branchEvent, branchId);

	// If node is extinct no edges, ancestor 1 edge, speciation 2 edges
	std::vector<PS::Edge*> edges = eventNodeTree->getEdgesToChildren();
	for(size_t iE=0; iE<edges.size(); ++iE) {
		itMapping_t itFind = edgeToDAG.find(edges[iE]->getId());
		assert(itFind != edgeToDAG.end());
		NodeDAG *prevNode = itFind->second;
		assert(prevNode != NULL);

		// Add an integration dag node for each branch leading to the speciation event
		/*NodeDAG *newNode = new NodeDAG(prevNode->getTime(), currentEvent->getTime(), edges[iE]->getId());
		prevNode->addParentNode(newNode);
		newNode->addChildNode(prevNode);*/
		NodeDAG *newNode = addSegmentedIntegrationNodes(prevNode->getTime(), currentEvent->getTime(), edges[iE]->getId(), prevNode);

		// This edge will collapse in the speciation, it's no longer active
		edgeToDAG[edges[iE]->getId()] = NULL;

		// Speciation node registering
		newNode->addParentNode(eventNodeDAG);
		eventNodeDAG->addChildNode(newNode);
	}

	// Add the speciation DAG node to the map
	if(eventNodeTree->isOriginNode()) {
		root = eventNodeDAG;
	} else {
		edgeToDAG[eventNodeTree->getEdgeToParent()->getId()] = eventNodeDAG;
		//std::cout << "Async edge : " << eventNodeTree->getEdgeToParent()->getId() << std::endl;
	}

	// If this node is an exinct tip -> then it's ready (but U is missing)
	if(eventNodeTree->isExtinct()) {
		tips.push_back(eventNodeDAG);
	}
}

void DAG::synchronousEvents(Event* currentEvent, mappingEdgeToDAG_t &edgeToDAG) {
	const std::list<PS::Edge*>&activeEdges = ptrScheduler->getActiveEdges(currentEvent->getTime());

	//std::cout << "Sync event : " << currentEvent->toString() << std::endl;

	// For each active edge
	for(std::list<PS::Edge*>::const_iterator it = activeEdges.begin(); it != activeEdges.end();  it++) {
		size_t idEdge = (*it)->getId();

		// Check if the edge identify an involved node
		std::vector<PS::Node*>::const_iterator itFindNode = std::find(currentEvent->getNodes().begin(), currentEvent->getNodes().end(), (*it)->getChild());
		if(itFindNode != currentEvent->getNodes().end()) continue;

		//std::cout << "Sync event [E] : " << (*it)->toString() << std::endl;
		//std::cout << "Sync event [C] : " << (*it)->getChild()->toString() << std::endl;
		//std::cout << "Sync event [P] : " << (*it)->getParent()->toString() << std::endl;

		// Get the active dag node
		itMapping_t itFind = edgeToDAG.find(idEdge);
		assert(itFind != edgeToDAG.end());
		NodeDAG *prevNode = itFind->second;
		assert(prevNode != NULL);

		// Create a new integration
		NodeDAG *newIntNode = addSegmentedIntegrationNodes(prevNode->getTime(), currentEvent->getTime(), (*it)->getId(), prevNode);
		/*NodeDAG *newIntNode = new NodeDAG(prevNode->getTime(), currentEvent->getTime(), (*it)->getId());
		prevNode->addParentNode(newIntNode);
		newIntNode->addChildNode(prevNode);*/

		// Create a new event node
		Event* branchEvent = new Event(currentEvent->getEvent(), currentEvent->getTime());
		NodeDAG *branchEventDAG = new NodeDAG(currentEvent->getTime(), branchEvent, (*it)->getId());
		newIntNode->addParentNode(branchEventDAG);
		branchEventDAG->addChildNode(newIntNode);

		// This edge will collapse in the speciation, it's no longer active
		edgeToDAG[idEdge] = branchEventDAG;
	}

	 if(currentEvent->checkEvent(SYNCHRONOUS_SAMPLING_EVENT) ||
		 currentEvent->checkEvent(SYCHRONOUS_SPECIATION_EVENT) ||
		 currentEvent->checkEvent(SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT)) {
		 // Deal with nodes
		 for(size_t iN=0; iN<currentEvent->getNodes().size(); ++iN) {
			PS::Node* eventNodeTree = currentEvent->getNodes()[iN];

			// Creating the dag node and mapping with previous DAG node if existing: add integration step too
			Event* branchEvent = new Event(currentEvent->getEvent(), eventNodeTree);
			NodeDAG *eventNodeDAG = new NodeDAG(currentEvent->getTime(), branchEvent, eventNodeTree->getEdgeToParent()->getId());

			// If node is extinct no edges, sampled ancesotr or speciation 2 edges
			std::vector<PS::Edge*> edges = eventNodeTree->getEdgesToChildren();
			assert(edges.empty() || edges.size() == 1 || edges.size() == 2);
			for(size_t iE=0; iE<edges.size(); ++iE) {
				NodeDAG *prevNode = edgeToDAG[edges[iE]->getId()];
				assert(prevNode != NULL);

				// Add an integration dag node for each branch leading to the speciation event
				NodeDAG *newNode = addSegmentedIntegrationNodes(prevNode->getTime(), currentEvent->getTime(), edges[iE]->getId(), prevNode);
				/*NodeDAG *newNode = new NodeDAG(prevNode->getTime(), currentEvent->getTime(), edges[iE]->getId());
				prevNode->addParentNode(newNode);
				newNode->addChildNode(prevNode);*/
				// This edge will collapse in the speciation, it's no longer active
				edgeToDAG[edges[iE]->getId()] = NULL;

				// New node registering
				newNode->addParentNode(eventNodeDAG);
				eventNodeDAG->addChildNode(newNode);
			}

			// Add the speciation DAG node to the map
			edgeToDAG[eventNodeTree->getEdgeToParent()->getId()] = eventNodeDAG;

			// If this node is an exinct tip -> then it's ready (but U is missing)
			if(eventNodeTree->isExtinct()) {
				tips.push_back(eventNodeDAG);
			}
		 }
	 }
}

void DAG::build() {

	clear();

	bool areEventPossible = true;
	for(size_t iE=0; iE<ptrScheduler->getEvents().size(); ++iE) {
		areEventPossible = areEventPossible && ptrScheduler->getEvents()[iE]->isEventPossible();
	}
	if(!areEventPossible) return;

	const std::vector<Event*> &events = ptrScheduler->getEvents();
	mappingEdgeToDAG_t edgeToDAG;

	assert(!events.empty());
	Event* currentEvent = events.front();

	// Init tips
	initTips(currentEvent, edgeToDAG);

	// For the next events:
	for(size_t iE=1; iE<events.size(); ++iE) {
		currentEvent = events[iE];

		if(currentEvent->checkEvent(NODE_EVENT) || currentEvent->checkEvent(FINAL_NODE_EVENT)) {
			asynchronousEvents(currentEvent, edgeToDAG);

		} else if(currentEvent->checkEvent(SYCHRONOUS_SPECIATION_EVENT) ||
					   currentEvent->checkEvent(SYNCHRONOUS_EXTINCTION_EVENT) ||
					   currentEvent->checkEvent(SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT) ||
					   currentEvent->checkEvent(SYNCHRONOUS_SAMPLING_EVENT) ||
					   currentEvent->checkEvent(SYNCHRONOUS_RATE_SHIFT) ||
					   currentEvent->checkEvent(SYNCHRONOUS_MONITORING_PROB)) {

			synchronousEvents(currentEvent, edgeToDAG);

		} else if(currentEvent->checkEvent(SYNCHRONOUS_RESCALING_EVENT)) {
			// Doing nothing here - we use a custom branchwise rescaling
		} else {
			assert(false && "Unknown event.");
		}
	}

	/*std::stringstream ss;
	recursiveToString(root, ss);
	std::cout << ss.str() << std::endl;*/
}

void DAG::rebuild() {
	clear();
	build();
}

void DAG::clear() {

	if(root != NULL) {
		root->recursiveClear();
		delete root;

		tips.clear();
		root = NULL;
	}
}

std::string DAG::toString() const {
	std::stringstream ss;
	recursiveToString(root, ss);
	return ss.str();
}

NodeDAG* DAG::getRoot() {
	return root;
}

const std::vector<NodeDAG*>& DAG::getTips() const {
	return tips;
}

void DAG::recursiveToString(const NodeDAG* aNode, std::stringstream &ss) const {
	ss << aNode->toString();
	ss << "------------------------------------------" << std::endl;

	for(size_t iN = 0; iN < aNode->getChildNodes().size(); ++iN) {
		recursiveToString(aNode->getChildNodes()[iN], ss);
	}
}

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */
