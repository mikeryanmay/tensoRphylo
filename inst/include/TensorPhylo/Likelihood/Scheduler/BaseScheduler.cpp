/*
 * BaseScheduler.cpp
 *
 *  Created on: Aug 23, 2019
 *      Author: xaviermeyer
 */

#include "BaseScheduler.h"

#include <limits>
#include <iterator>
#include <iostream>
#include <vector>

#include "Event.h"
#include "Data/Structure/IncTreeStructure.h"
#include "Parameters/IncParameterContainer.h"
#include "SynchronousEvents/IncSynchronousEvents.h"

namespace Likelihood {
namespace Scheduler {

const double BaseScheduler::MAX_SEGMENT_SIZE_WITHOUT_RESCALING = 0.8;

BaseScheduler::BaseScheduler(PS::TreeSharedPtr aPtrTree) : ptrTree(aPtrTree) {
	initEvents();
	updated = true;
}

BaseScheduler::~BaseScheduler() {
	clearEvents();
}


void BaseScheduler::setMonitoringProbes(const std::vector<double> &aTimes) {
	for(size_t iT=0; iT<aTimes.size(); ++iT) {
		assert(aTimes[iT] < events.back()->getTime() && "A monitoring probe is defined after the origination age.");
		addSynchronousEvent(aTimes[iT], SYNCHRONOUS_MONITORING_PROB);
	}
	updated = true;
}

void BaseScheduler::setRateShiftEvents(Parameters::AsyncContainerSharedPtr ptrAsyncParametersContainer) {

	addRateShiftEvents(ptrAsyncParametersContainer->timesLambda, SYNCHRONOUS_RATE_SHIFT);
	addRateShiftEvents(ptrAsyncParametersContainer->timesMu, SYNCHRONOUS_RATE_SHIFT);
	addRateShiftEvents(ptrAsyncParametersContainer->timesDelta, SYNCHRONOUS_RATE_SHIFT);
	addRateShiftEvents(ptrAsyncParametersContainer->timesPhi, SYNCHRONOUS_RATE_SHIFT);
	addRateShiftEvents(ptrAsyncParametersContainer->timesEta, SYNCHRONOUS_RATE_SHIFT);
	addRateShiftEvents(ptrAsyncParametersContainer->timesOmega, SYNCHRONOUS_RATE_SHIFT);
	updated = true;
}

void BaseScheduler::setSynchronousEvents(SynchronousEvents::ContainerSharedPtr ptrSyncEventContainer) {

	{ // Speciation
		std::vector<double> &times = ptrSyncEventContainer->getPtrMassSpeciation()->getEventTimes();
		for(size_t iT=0; iT<times.size(); ++iT) {
			addMassSpeciationEvent(times[iT]);
		}
	}

	{ // Extinction
		std::vector<double> &times = ptrSyncEventContainer->getPtrMassExtinction()->getEventTimes();
		for(size_t iT=0; iT<times.size(); ++iT) {
			addMassExtinctionEvent(times[iT]);
		}
	}

	{ // Sampling
		std::vector<double> &times = ptrSyncEventContainer->getPtrMassSampling()->getEventTimes();
		assert(times.size() >= 1 && "The initial (t=0) mass sampling event must be defined as a mass sampling event.");
		assert(times[0] == 0. && "The initial (t=0) mass sampling event must be defined and provide the initial sampling probabilities (Rho_ij where i is state and j ).");
		for(size_t iT=1; iT<times.size(); ++iT) { // We skip the initial mass sampling event (at time 0) as it is the initial condition
			addMassSamplingEvent(times[iT]);
		}
	}

	{ // Destructive sampling
		std::vector<double> &times = ptrSyncEventContainer->getPtrMassDestrSampling()->getEventTimes();
		for(size_t iT=0; iT<times.size(); ++iT) {
			addMassDestructiveSamplingEvent(times[iT]);
		}
	}
	updated = true;
}

void BaseScheduler::addMassExtinctionEvent(double aTime) {
	addSynchronousEvent(aTime, SYNCHRONOUS_EXTINCTION_EVENT);
	updated = true;
}

void BaseScheduler::addMassSamplingEvent(double aTime) {
	addSynchronousEvent(aTime, SYNCHRONOUS_SAMPLING_EVENT);
	updated = true;
}

void BaseScheduler::addMassSpeciationEvent(double aTime) {
	addSynchronousEvent(aTime, SYCHRONOUS_SPECIATION_EVENT);
	updated = true;
}


void BaseScheduler::addMassDestructiveSamplingEvent(double aTime) {
	addSynchronousEvent(aTime, SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT);
	updated = true;
}

void BaseScheduler::removeRateShiftEvents() {
	removeEventsByType(SYNCHRONOUS_RATE_SHIFT);
	updated = true;
}

void BaseScheduler::removeSynchronousEvents() {
	removeEventsByType(SYNCHRONOUS_EXTINCTION_EVENT);
	removeEventsByType(SYNCHRONOUS_SAMPLING_EVENT);
	removeEventsByType(SYCHRONOUS_SPECIATION_EVENT);
	removeEventsByType(SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT);
	updated = true;
}

void BaseScheduler::removeRescalingEvents() {
	removeEventsByType(SYNCHRONOUS_RESCALING_EVENT);
	updated = true;
}

void BaseScheduler::removeEventsByType(eventType_t aEventType) {

	std::vector<Event*>::iterator it = events.begin();

	while(it != events.end()) { // Looping through event and removing async. rate shifts
		if((*it)->checkEvent(aEventType)) {

			// Finding the position of the event
			size_t pos = std::distance(events.begin(), it);
			// Erase the associated layer of edges
			layeredEdges.erase(layeredEdges.begin() + (pos - 1));
			// Erase the event
			delete (*it);
			it = events.erase(it);

		} else {
			it++;
		}
	}

}

const std::vector<Event*>& BaseScheduler::getEvents() const {
	return events;
}

const std::list<PS::Edge*>& BaseScheduler::getActiveEdges(double time) const {


	assert(!events.empty() && "[getActiveEdges] : No event registered.");
	assert(time < events.back()->getTime() && "[getActiveEdges] : This function is only defined from 0 to age of root.");

	for(size_t iE=1; iE<events.size(); ++iE) {
		if(time < events[iE]->getTime()) {
			return layeredEdges[iE-1];
		}
	}

	assert(false && "[getActiveEdges] : No layer of edges found but time is smaller than age of root.");
	return layeredEdges.back();
}

PS::TreeSharedPtr BaseScheduler::getPtrTree() const {
	return ptrTree;
}

bool BaseScheduler::hasBeenUpdated() const {
	return updated;
}

void BaseScheduler::clearHasBeenUpdatedFlag() {
	updated = false;
}

void BaseScheduler::initEvents() {

	if(!events.empty()) {
		clearEvents();
	}

	 // Starting point: leaf nodes
	const std::vector<PS::Node*> &extantNodes = ptrTree->getExtantNodes();
	Event* startEvent = new Event(PRESENT_TIME_EVENT, extantNodes);
	events.push_back(startEvent);

	while(!events.back()->checkEvent(FINAL_NODE_EVENT)) {

		// Insert the next layer of edges and return the next event (internal / root nodes)
		std::vector< PS::Node* > nextEventNodes = defineNextEdgesLayerAndEvent(events.back());
		assert(!nextEventNodes.empty());

		Event* nextEvent = NULL;
		if(nextEventNodes.size() == 1 && nextEventNodes.front()->isOriginNode()) {
			nextEvent = new Event(FINAL_NODE_EVENT, nextEventNodes.front());
		} else { // Ancestral / fossil
			// Sanity check FIXME remove when sure everything is ok
			for(size_t iN=0; iN<nextEventNodes.size(); ++iN) {
				assert(!nextEventNodes[iN]->isOriginNode());
			}

			nextEvent = new Event(NODE_EVENT, nextEventNodes);
		}
		events.push_back(nextEvent);
	}
}

void BaseScheduler::defineAndSetRescalingEvents() {

	removeRescalingEvents();

	std::vector< std::pair<double, size_t> > rescalingTime; // Trick to get the argsort
	std::vector<Phylogeny::Structure::Node*> nodesId;
	std::vector<bool> associatedToOtherEvent;

	// For each edges
	const std::vector<Phylogeny::Structure::Edge*>& edges = ptrTree->getEdges();
	for(size_t iE=0; iE<edges.size(); ++iE) {
		// If the edges is longer than a given length, we request for rescaling events
		if(edges[iE]->getLength() > MAX_SEGMENT_SIZE_WITHOUT_RESCALING) {
			double edgeLength = edges[iE]->getLength();
			size_t nSubSegment = std::floor(edgeLength/MAX_SEGMENT_SIZE_WITHOUT_RESCALING);
			double subsegmentSize = edgeLength/nSubSegment;

			//std::cout << "Treating edge : " << edges[iE]->toString();
			//std::cout << "With child node : " << edges[iE]->getChild()->toString();

			// Request a rescaling at time child_node + i*subsegment_length
			for(size_t iS=1; iS<nSubSegment; ++iS) {
				double requestedTime = edges[iE]->getChild()->getAge()+subsegmentSize*iS;
				//std::cout << "Subsegment boundary " << iS << " with requested time = " << requestedTime << std::endl;

				int firstEventAfterReqTime = getFirstEventAfterTimeT(requestedTime);
				assert(firstEventAfterReqTime >= 0);

				int firstEventPriorReqTime = firstEventAfterReqTime - 1;
				assert(firstEventAfterReqTime >= 0);

				// if it is within a 30% subsegment length tolerance of the requested time, we update the req time
				if(fabs(events[firstEventAfterReqTime]->getTime() - requestedTime) < (0.3*MAX_SEGMENT_SIZE_WITHOUT_RESCALING)/2.) {
					requestedTime = events[firstEventAfterReqTime]->getTime();
					associatedToOtherEvent.push_back(true);
					//std::cout << "Updated to next event  :" << events[firstEventAfterReqTime]->toString() << std::endl;
				} else if(fabs(events[firstEventPriorReqTime]->getTime() - requestedTime) < (0.3*MAX_SEGMENT_SIZE_WITHOUT_RESCALING)/2.) {
					requestedTime = events[firstEventPriorReqTime]->getTime();
					//std::cout << "Updated to prior event  :" << events[firstEventPriorReqTime]->toString() << std::endl;
					associatedToOtherEvent.push_back(true);
				} else {
					associatedToOtherEvent.push_back(false);
				}

				nodesId.push_back(edges[iE]->getChild());
				rescalingTime.push_back(std::make_pair(requestedTime, rescalingTime.size()));
				//std::cout << "Final time  :" << requestedTime << std::endl;
			}
			//std::cout << "--------------------------------------------" << std::endl;
		}
	}

	if(rescalingTime.empty()) return; // no rescaling times, we are out

	updated = true;

	// Otherwise, regroup them and add events
	// Order and regroup events
	std::sort(rescalingTime.begin(), rescalingTime.end()); // Getting the sorted time and argsort

	// First element
	double firstTime = -1;
	std::vector<Phylogeny::Structure::Node*> rescalingEventNodes;
	for(size_t iR=0; iR<rescalingTime.size(); ++iR) {

		double time = rescalingTime[iR].first;
		if(firstTime < 0.) firstTime = time;

		size_t pos = rescalingTime[iR].second;

		rescalingEventNodes.push_back(nodesId[pos]);

		// There is 3 reasons to create an event at this point:
		// 1) We are the last rescaling event
		// 2) There is already a synchronous event at this time and the next rescaling event differ in time
		// 3) The difference between rescaling events unrelated to an existing synchronous event exceed 30% of the MAX_SEGMENT_SIZE
		if(iR == rescalingTime.size()-1 ||
		   (associatedToOtherEvent[pos] && rescalingTime[iR+1].first != time) ||
				(rescalingTime[iR+1].first - firstTime) > 0.3*MAX_SEGMENT_SIZE_WITHOUT_RESCALING) {

			if(!associatedToOtherEvent[pos]) { // if we are not associated, we take the mean time
				time = 0.;
				for(size_t iS=0; iS<rescalingEventNodes.size(); ++iS) {
					time += rescalingTime[iR-iS].first/rescalingEventNodes.size();
				}
				//std::cout << "We are a bunch of free rescaling, thus averaging time to t = " << time << std::endl;
			} else {
				//std::cout << "We are associated with an event, thus keeping time t = " << time << std::endl;
			}

			// We create and insert an event
			int index = getFirstEventAfterTimeT(time);
			Event* newEvent = new Event(SYNCHRONOUS_RESCALING_EVENT, time, rescalingEventNodes);
			//std::cout << "Creating event : " << newEvent->toString() << std::endl;
			std::vector<Event*>::iterator itE = events.begin();
			std::advance(itE, index);
			events.insert(itE, newEvent);

			// update the edges layer
			std::vector< edgesList_t >::iterator itL = layeredEdges.begin();
			std::advance(itL, index-1);
			edgesList_t layerCopy  = *itL;
			layeredEdges.insert(itL, layerCopy);

			// We clear the rescalingEventNode buffer
			rescalingEventNodes.clear();
			firstTime = -1.;
		}
	}
	//std::cout << "--------------------------------------------" << std::endl;
}


void BaseScheduler::clearEvents() {
	for(size_t iE=0; iE<events.size(); ++iE) {
		delete events[iE];
	}
	events.clear();
}

std::vector<PS::Node*> BaseScheduler::defineNextEdgesLayerAndEvent(Event *lastEvent) {

	const std::vector<PS::Node*> &currentNodes = lastEvent->getNodes();

	// Find the next layer of edge
	edgesList_t nextEdges;
	if(lastEvent->checkEvent(PRESENT_TIME_EVENT)) {
		// Current nodes contains all the leaf nodes
		// 1) we must add all edges
		// 2) we must find the next candidate ancestral node (event)
		for(size_t iN=0; iN<currentNodes.size(); ++iN) {
			PS::Edge *edge = currentNodes[iN]->getEdgeToParent();

			assert(edge);
			nextEdges.push_back(edge);
		}
	} else if(lastEvent->checkEvent(NODE_EVENT)) {
		// Current nodes contains the ancestral node(s)

		// Starting from last edges
		nextEdges = layeredEdges.back();

		// For each ancestral node with the same age
		for(size_t iN=0; iN<currentNodes.size(); ++iN){
			// 1) Remove old child edges (only for ancestral, not for fossil)
			std::vector<PS::Edge*> childrenEdges = currentNodes[iN]->getEdgesToChildren();
			for(size_t iC=0; iC<childrenEdges.size(); ++iC) {
				nextEdges.remove(childrenEdges[iC]);
			}
			// 2) Add new parent edge
			nextEdges.push_back(currentNodes[iN]->getEdgeToParent());
		}
	} else if(lastEvent->checkEvent(FINAL_NODE_EVENT)) {
		// Current nodes contains the root node
		// Check that we only have the last two edges as the last layer
		assert(layeredEdges.back().size() == 2);
		assert(false); // We should not arrive here

		// force an empty return
		std::vector<PS::Node*> empty;
		return empty;
	}

	// Add the new edges to the layers
	layeredEdges.push_back(nextEdges);

	// Find the next event(s) :
	// 1) next speciation event
	std::vector<PS::Node*> nextEventNodes;
	double nextEventAge = std::numeric_limits<double>::max();
	for(itEdgesList_t itE=nextEdges.begin(); itE != nextEdges.end(); ++itE) {
		PS::Edge *edge = (*itE);

		if(edge->getParent()->getAge() < nextEventAge) { // If sooner that current event, clear and memorize
			nextEventAge = edge->getParent()->getAge();
			nextEventNodes.clear();
			nextEventNodes.push_back(edge->getParent());
		} else if(edge->getParent()->getAge() == nextEventAge && // If same age and different node, add
					  std::find(nextEventNodes.begin(), nextEventNodes.end(), edge->getParent()) == nextEventNodes.end() ) {
			nextEventNodes.push_back(edge->getParent());
		}
	}

	//2 next fossil event
	double currentEventAge = events.back()->getTime();
	const std::vector<PS::Node*>& fossilNodes = ptrTree->getExtinctNodes();
	for(size_t iN=0; iN<fossilNodes.size(); ++iN) {
		if(fossilNodes[iN]->getAge() <= currentEventAge) continue; // Already dealt with

		if(fossilNodes[iN]->getAge() < nextEventAge) { // If sooner that current event, clear and memorize
			nextEventAge = fossilNodes[iN]->getAge();
			nextEventNodes.clear();
			nextEventNodes.push_back(fossilNodes[iN]);
		} else if(fossilNodes[iN]->getAge() == nextEventAge &&  // If same age and different node, add
				      std::find(nextEventNodes.begin(), nextEventNodes.end(), fossilNodes[iN]) == nextEventNodes.end() ) {
			nextEventNodes.push_back(fossilNodes[iN]);
		}
	}

	return nextEventNodes;
}


int BaseScheduler::getFirstEventAfterTimeT(double aTime) {
	for(size_t iE=0; iE < events.size(); ++iE) {
		//std::cout << "[GET] Event time = " << aTime << " vs event time =  " << events[iE]->getTime() << " -- " << std::scientific << fabs(aTime - events[iE]->getTime()) << (aTime == events[iE]->getTime() ? " -- SAME!" : "-- not same??") << std::endl;
		if(aTime <= events[iE]->getTime()) {
			return iE;
		} else if (std::fabs(aTime - events[iE]->getTime()) <= 2.*std::numeric_limits<double>::epsilon()) {
			return iE;
		}
	}
	return -1;
}


void BaseScheduler::addSynchronousEvent(double aTime, eventType_t aEventType) {

	int index = getFirstEventAfterTimeT(aTime);

	//std::cout << "Found : " << events[index]->toString() << std::endl;
	//std::cout << "Event time = " << aTime << " vs event time =  " << events[index]->getTime() << (aTime == events[index]->getTime() ? " -- SAME!" : "-- not same??") << std::endl;

	bool areSameUpToEpsilon = (std::fabs(aTime - events[index]->getTime()) <= 2.*std::numeric_limits<double>::epsilon());
	if(areSameUpToEpsilon && aEventType != SYNCHRONOUS_MONITORING_PROB) {
		//std::cout << "Updating " << std::endl;
		events[index]->updateEvent(aEventType);
	} else { // No event at this time, inserting a new one.
		//std::cout << "New " << std::endl;
		Event* newEvent = new Event(aEventType, aTime);
		std::vector<Event*>::iterator itE = events.begin();
		std::advance(itE, index);
		events.insert(itE, newEvent);

		std::vector< edgesList_t >::iterator itL = layeredEdges.begin();
		std::advance(itL, index-1);
		edgesList_t layerCopy  = *itL;
		layeredEdges.insert(itL, layerCopy);
	}
}

void BaseScheduler::addRateShiftEvents(const std::vector<double> &aTimes, eventType_t aEventType) {
	assert(aEventType == SYNCHRONOUS_RATE_SHIFT);

	// We want to insure that there is an event at the RATE_SHIFT time(s)
	for(size_t iT=0; iT < aTimes.size(); ++iT) {

		double time = aTimes[iT];
		int index = getFirstEventAfterTimeT(time);

		if(events[index]->getTime() == time) { // There is already an event
			return;
		} else { // No event at this time, inserting a new one.
			Event* newEvent = new Event(aEventType, time);
			std::vector<Event*>::iterator itE = events.begin();
			std::advance(itE, index);
			events.insert(itE, newEvent);

			std::vector< edgesList_t >::iterator itL = layeredEdges.begin();
			std::advance(itL, index-1);
			edgesList_t layerCopy  = *itL;
			layeredEdges.insert(itL, layerCopy);
		}
	}
}

} /* namespace Scheduler */
} /* namespace Likelihood */
