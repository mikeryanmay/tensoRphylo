/*
 * Event.cpp
 *
 *  Created on: Aug 23, 2019
 *      Author: xaviermeyer
 */

#include "Event.h"

#include <cassert>
#include <sstream>
#include <iostream>
#include <boost/assign/list_of.hpp>

#include "Data/Structure/Node.h"
#include "Data/Structure/Edge.h"

namespace Likelihood {
namespace Scheduler {

Event::Event(eventType_t aEType, double aTime) : eventType(aEType), time(aTime) {
	isImpossible = false;
}

Event::Event(eventType_t aEType, PS::Node *aNode) : eventType(aEType), eventNodes(1, aNode) {
	isImpossible = false;
	assert(aNode);
	time = eventNodes.front()->getAge();
}

Event::Event(eventType_t aEType, const std::vector<PS::Node*> &aNodes) : eventType(aEType), eventNodes(aNodes) {
	isImpossible = false;
	assert(!eventNodes.empty() && eventNodes.front());
	time = eventNodes.front()->getAge();
	for(size_t iN=0; iN<eventNodes.size(); ++iN) {
		assert(eventNodes[iN]);
		assert(eventNodes[iN]->getAge() == time);
	}
}

Event::Event(eventType_t aEType, double aTime, const std::vector<PS::Node*> &aNodes) : eventType(aEType), time(aTime), eventNodes(aNodes) {
	isImpossible = false;
	assert(!eventNodes.empty());
	assert(aEType == SYNCHRONOUS_RESCALING_EVENT);
	for(size_t iN=0; iN<eventNodes.size(); ++iN) {
		assert(eventNodes[iN]);
		assert(!eventNodes[iN]->isOriginNode());
		assert(eventNodes[iN]->getAge() < aTime);
		assert(eventNodes[iN]->getParentNode()->getAge() > aTime);
	}
}

Event::~Event() {
}

void Event::updateEvent(eventType_t aEventType) {
	// FIXME what happens if we update the root event ? Is that even possible ?
	std::vector<eventType_t> synchronousEvent = boost::assign::list_of(SYCHRONOUS_SPECIATION_EVENT)(SYNCHRONOUS_EXTINCTION_EVENT)(SYNCHRONOUS_SAMPLING_EVENT)(SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT);

	bool isOldSynchronousEvent = std::find(synchronousEvent.begin(), synchronousEvent.end(), eventType) != synchronousEvent.end();
	bool isNewSynchronousEvent = std::find(synchronousEvent.begin(), synchronousEvent.end(), aEventType) != synchronousEvent.end();
	if(isOldSynchronousEvent && isNewSynchronousEvent) {
		isImpossible = true;
	}
	eventType = aEventType;
}

bool Event::checkEvent(eventType_t aEventType) const {
	return eventType == aEventType;
}

eventType_t Event::getEvent() const {
	return eventType;
}

void Event::setTime(double aTime) {
	time = aTime;
}

double Event::getTime() const {
	return time;
}

const std::vector<PS::Node*>& Event::getNodes() const {
	return eventNodes;
}


bool Event::isEventPossible() const {

	// First we check if two synchronous event have been defined at the same time
	if(isImpossible) return false;

	// The depending on the event we have to check the nodes registered
	bool isValid = true;
	switch (eventType) {
		case PRESENT_TIME_EVENT:
			isValid = !eventNodes.empty(); 						// We must have a least one extant node
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {				// ... and all nodes must be extant
				isValid = isValid && eventNodes[iN]->isExtant();
			}
			break;
		case NODE_EVENT:
				// We must have one node of type speciation, sampled ancestor or speciation
			isValid = eventNodes.size() == 1 &&
						  (eventNodes[0]->isExtinct() ||
						   eventNodes[0]->isSampledAncestor() ||
						   eventNodes[0]->isSpeciationNode());
			if ( isValid == false ) {
				std::cout << "this speciation node has " << eventNodes.size() << " nodes." << std::endl;
				std::cout << "this speciation node is at time " << getTime() << std::endl;
			}
			break;
		case SYCHRONOUS_SPECIATION_EVENT:
			// Nodes must be speciation nodes:
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				isValid = isValid && eventNodes[iN]->isSpeciationNode();
			}
			break;
		case SYNCHRONOUS_EXTINCTION_EVENT:
			isValid = eventNodes.empty();							// There must be no nodes
			break;
		case SYNCHRONOUS_SAMPLING_EVENT:
			// Nodes must be either extinct node or sampled ancestor
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				isValid = isValid &&
						(eventNodes[iN]->isExtinct() ||
						eventNodes[iN]->isSampledAncestor());
			}
			break;
		case SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT:
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {				// Nodes must be :
				isValid = isValid &&
						eventNodes[iN]->isExtinct();				 				// extinct node
			}
			break;
		case FINAL_NODE_EVENT:													// Final node must be
			isValid = eventNodes.size() == 1 &&
			eventNodes.front()->isOriginNode();								// unique and the origination node (which can be the root node)
			break;
		case SYNCHRONOUS_RATE_SHIFT:
			break;
		case SYNCHRONOUS_MONITORING_PROB:
			break;
		case SYNCHRONOUS_RESCALING_EVENT:
			break;
		default:
			assert(false && "Unregistered event.");
			break;
	}

	if ( isValid == false ) {
		std::cout << "This " << eventType << " is invalid." << std::endl;
	}

	return isValid;
}

std::string Event::toString() const {

	std::stringstream ss;

	switch (eventType) {
		case PRESENT_TIME_EVENT:
			ss << "Present time (all extant taxa).";
			ss << std::endl << "Nodes : [ ";
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				ss << eventNodes[iN]->getId() << " ,";
			}
			ss << "]";
			break;
		case NODE_EVENT:
			assert(!eventNodes.empty());
			ss << "Internal node (fossil, internal, root if origin != root) : [ ";
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				ss << eventNodes[iN]->getId() << " ,";
			}
			ss << "] - age = " << time;
			break;
		case SYCHRONOUS_SPECIATION_EVENT:
			ss << "Synchronous speciation event : " << "  time = " << time;
			ss << std::endl << "Nodes : [ ";
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				ss << eventNodes[iN]->getId() << " ,";
			}
			ss << "]";
			break;
		case SYNCHRONOUS_EXTINCTION_EVENT:
			assert(eventNodes.empty());
			ss << "Synchronous extinction event : " << "  time = " << time;
			break;
		case SYNCHRONOUS_SAMPLING_EVENT:
			ss << "Synchronous sampling event : " << "  time = " << time;
			ss << std::endl << "Nodes : [ ";
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				ss << eventNodes[iN]->getId() << " ,";
			}
			ss << "]";
			break;
		case SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT:
			ss << "Synchronous destructive sampling event : " << "  time = " << time;
			ss << std::endl << "Nodes : [ ";
			for(size_t iN=0; iN<eventNodes.size(); ++iN) {
				ss << eventNodes[iN]->getId() << " ,";
			}
			ss << "]";
			break;
		case FINAL_NODE_EVENT:
			assert(!eventNodes.empty() && eventNodes.size() == 1);
			ss << "Final node event (origin/root) : " << eventNodes.front()->getId() << " - age = " << time;
			break;
		case SYNCHRONOUS_RATE_SHIFT:
			ss << "Asynchronous rate shift at time t = " << time;
			break;
		case SYNCHRONOUS_MONITORING_PROB:
			ss << "Monitoring probe at time t = " << time;
			break;
		case SYNCHRONOUS_RESCALING_EVENT:
			ss << "Rescaling at time t=" << time << " for branches: ";
			for(size_t iN=0; iN<eventNodes.size(); ++iN) ss << eventNodes[iN]->getEdgeToParent()->getId() << ", ";
			break;
		default:
			assert(false && "Unregistered event.");
			break;
	}

	ss << (isEventPossible() ? " \t possible event" : "\t impossible event");

	return ss.str();
}

} /* namespace Scheduler */
} /* namespace Likelihood */
