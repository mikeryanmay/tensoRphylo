/*
 * NodeDAG.cpp
 *
 *  Created on: Apr 13, 2020
 *      Author: xaviermeyer
 */

#include "NodeDAG.h"

#include "../Event.h"
#include "Data/Structure/Node.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>

namespace Likelihood {
namespace Scheduler {
namespace DAG {


size_t NodeDAG::idSeq = 0;

NodeDAG::NodeDAG(double aTime, Event *aEvent, long int aIdEdge) :
		id(idSeq++), taskType(EVENT_TASK), event(aEvent),
		 idIntEdge(aIdEdge), startTime(aTime), endTime(aTime) {

}

NodeDAG::NodeDAG(double aStartTime, double aEndTime, long int aIdEdge):
				id(idSeq++), taskType(INTEGRATE_TASK), event(NULL),
				idIntEdge(aIdEdge), startTime(aStartTime), endTime(aEndTime) {

}

NodeDAG::~NodeDAG() {

	if(event != NULL) {
		delete event;
		event = NULL;
	}

	assert(probStates.empty());
}

void NodeDAG::attachState(StateType *aPState) {
	probStates.push_back(aPState);
}

void NodeDAG::detachState(StateType *aPState) {
	itVecStates_t itFind = std::find(probStates.begin(), probStates.end(), aPState);
	assert(itFind != probStates.end() && "Detaching an unattached state.");

	probStates.erase(itFind);
}

void NodeDAG::detachStates() {
	probStates.clear();
}


std::vector<StateType*>& NodeDAG::getProbStates() {
	return probStates;
}

const std::vector<StateType*>& NodeDAG::getProbStates() const {
	return probStates;
}

bool NodeDAG::isReady() const {

	bool ready = true;
	for(size_t iC=0; iC<childNodes.size(); ++iC) {
		ready = ready && childNodesReady[iC];
	}
	return ready;

}

bool NodeDAG::isIntegrationTask() const {
	return taskType == INTEGRATE_TASK;
}

bool NodeDAG::isEventTask() const {
	return taskType == EVENT_TASK;
}


size_t NodeDAG::getId() const {
	return id;
}

double NodeDAG::getStartTime() const {
	assert(taskType != EVENT_TASK && startTime != endTime);
	return startTime;
}

double NodeDAG::getEndTime() const {
	assert(taskType == INTEGRATE_TASK && startTime != endTime);
	return endTime;
}

double NodeDAG::getTime() const {
	assert(taskType == EVENT_TASK && startTime == endTime);
	return startTime;
}


Event* NodeDAG::getEvent() {
	assert(taskType == EVENT_TASK && event != NULL);
	return event;
}

void NodeDAG::addChildNode(NodeDAG* aNode) {
	childNodes.push_back(aNode);
	childNodesReady.push_back(false);
}

void NodeDAG::addParentNode(NodeDAG* aNode) {
	parentNodes.push_back(aNode);
}

void NodeDAG::removeChildNode(NodeDAG* aNode) {
	std::vector<NodeDAG*>::iterator itFind = std::find(childNodes.begin(), childNodes.end(), aNode);
	assert(itFind != childNodes.end() && "Removing an unregistered dag node.");

	childNodesReady.erase(childNodesReady.begin() + std::distance(childNodes.begin(), itFind));
	childNodes.erase(itFind);
}

void NodeDAG::removeParentNode(NodeDAG* aNode) {
	std::vector<NodeDAG*>::iterator itFind = std::find(parentNodes.begin(), parentNodes.end(), aNode);
	assert(itFind != parentNodes.end() && "Removing an unregistered dag node.");

	parentNodes.erase(itFind);
}

long int NodeDAG::getEdgeIdMapping() const {
	return idIntEdge;
}

void NodeDAG::setChildNodeReady(NodeDAG* aChildNode) {
	std::vector<NodeDAG*>::iterator itFind = std::find(childNodes.begin(), childNodes.end(), aChildNode);
	assert(itFind != childNodes.end() && "Removing an unregistered dag node.");

	childNodesReady[std::distance(childNodes.begin(), itFind)] = true;
}

void NodeDAG::signalReadyToParents() {
	// Signal to parent
	for(size_t iP=0; iP<parentNodes.size(); ++iP) {
		parentNodes[iP]->setChildNodeReady(this);
	}
}

void NodeDAG::resetChildReadyState() {
	// Reset state
	bool wasIReady = true;
	for(size_t iC=0; iC<childNodes.size(); ++iC) {
		wasIReady = wasIReady && childNodesReady[iC];
		childNodesReady[iC] = false;
	}
	assert(wasIReady);
}

std::vector<NodeDAG*>& NodeDAG::getChildNodes() {
	return childNodes;
}

const std::vector<NodeDAG*>& NodeDAG::getChildNodes() const {
	return childNodes;
}

std::vector<NodeDAG*>& NodeDAG::getParentNodes() {
	return parentNodes;
}

const std::vector<NodeDAG*>& NodeDAG::getParentNodes() const {
	return parentNodes;
}


std::string NodeDAG::toString() const {
	std::stringstream ss;
	ss << "ID : " << id << " -- Event: " << std::endl;;
	if(isEventTask()) {
		ss << event->toString() << std::endl;
	} else {
		ss << " integration from " << startTime << " to " << endTime << " ( edge " << idIntEdge << " )" << std::endl;
	}
	ss << "Parent nodes : ";
	for(size_t iN  = 0; iN < parentNodes.size(); ++iN) {
		ss << parentNodes[iN]->getId() << "\t";
	}
	ss << std::endl;
	ss << "Child nodes : ";
	for(size_t iN  = 0; iN < childNodes.size(); ++iN) {
		ss << childNodes[iN]->getId() << (childNodesReady[iN] ? "[R]": "[!R]") << "\t";
	}
	ss << std::endl;

	return ss.str();
}


void NodeDAG::identifyTips(std::vector<NodeDAG*> &tips, std::set<size_t> &visitedTips) {

	if(childNodes.empty()) { // If we are a tip
		if(visitedTips.count(id) == 0) { // and have not been visited by anybody
			tips.push_back(this);
		}
		visitedTips.insert(id);
	} else {
		for(size_t iC=0; iC < childNodes.size(); ++iC) {
			childNodes[iC]->identifyTips(tips, visitedTips);
		}
	}

}


double NodeDAG::identifySubDAGSumBL() {
	return identifySubDAGSumBL(this);
}

double NodeDAG::identifySubDAGSumBL(NodeDAG* aNode) {
	assert(aNode);
	double sumBL = aNode->endTime - aNode->startTime;
	for(size_t i=0; i<aNode->getChildNodes().size(); ++i) {
		sumBL += identifySubDAGSumBL(aNode->getChildNodes()[i]);
	}
	return sumBL;
}

void NodeDAG::recursiveClear() {

	// Delete child
	std::vector<NodeDAG*> tmpChildNodes = childNodes;
	for(size_t iC=0; iC<tmpChildNodes.size(); ++iC) {
		tmpChildNodes[iC]->recursiveClear();
		delete tmpChildNodes[iC];
	}

	// Unregister from parent (required to avoid double delete since we are a DAG and not a tree)
	for(size_t iP=0; iP<parentNodes.size(); ++iP) {
		parentNodes[iP]->removeChildNode(this);
	}
}

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */
