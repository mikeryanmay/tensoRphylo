/*
 * BaseScheduler.h
 *
 *  Created on: Aug 23, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_BASESCHEDULER_H_
#define LIKELIHOOD_SCHEDULER_BASESCHEDULER_H_

#include <boost/smart_ptr/shared_ptr.hpp>
#include <list>
#include <vector>

#include "EventType.h"
#include "Data/Structure/IncFwdTreeStructure.h"
#include "Parameters/IncFwdParameterContainer.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"

namespace Likelihood {
namespace Scheduler {

class Event;

class BaseScheduler {
public:
	static const double MAX_SEGMENT_SIZE_WITHOUT_RESCALING;
public:
	BaseScheduler(PS::TreeSharedPtr aPtrTree);
	virtual ~BaseScheduler();

	void setMonitoringProbes(const std::vector<double> &aTimes);
	void setRateShiftEvents(Parameters::AsyncContainerSharedPtr ptrAsyncParametersContainer);
	void setSynchronousEvents(SynchronousEvents::ContainerSharedPtr ptrSyncEventContainer);
	void defineAndSetRescalingEvents();

	void addMassExtinctionEvent(double aTime);
	void addMassSamplingEvent(double aTime);
	void addMassSpeciationEvent(double aTime);
	void addMassDestructiveSamplingEvent(double aTime);

	void removeRateShiftEvents();
	void removeSynchronousEvents();
	void removeRescalingEvents();

	const std::vector<Event*>& getEvents() const;

	// If argument time is the same as an event, the edges layer for t+dt is returned.
	const std::list<PS::Edge*>& getActiveEdges(double time) const;

	PS::TreeSharedPtr getPtrTree() const;

	bool hasBeenUpdated() const;
	void clearHasBeenUpdatedFlag();

protected:

	bool updated;

	PS::TreeSharedPtr ptrTree;

	std::vector<Event*> events;

	typedef std::list<PS::Edge*> edgesList_t;
	typedef edgesList_t::iterator itEdgesList_t;
	std::vector< edgesList_t > layeredEdges;

	void initEvents();
	std::vector<PS::Node*> defineNextEdgesLayerAndEvent(Event *lastEvent);
	void clearEvents();

	int getFirstEventAfterTimeT(double aTime);

	void addSynchronousEvent(double aTime, eventType_t aEventType);
	void addRateShiftEvents(const std::vector<double> &aTimes, eventType_t aEventType);

	void removeEventsByType(eventType_t aEventType);
};

} /* namespace Scheduler */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SCHEDULER_BASESCHEDULER_H_ */
