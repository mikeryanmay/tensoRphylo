/*
 * Event.h
 *
 *  Created on: Aug 23, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_EVENT_H_
#define LIKELIHOOD_SCHEDULER_EVENT_H_

#include <string>
#include <vector>

#include "EventType.h"

namespace Phylogeny {
namespace Structure {
	class Node;
}
}

namespace Likelihood {
namespace Scheduler {

namespace PS = ::Phylogeny::Structure;

class Event {

public:
	Event(eventType_t aEType, double aTtime);
	Event(eventType_t aEType, PS::Node *aNode);

	Event(eventType_t aEType, const std::vector<PS::Node*> &aNodes);
	Event(eventType_t aEType, double aTime, const std::vector<PS::Node*> &aNodes);
	virtual ~Event();

	void updateEvent(eventType_t aEventType);
	bool checkEvent(eventType_t aEventType) const;
	eventType_t getEvent() const;

	void setTime(double aTime);
	double getTime() const;
	const std::vector<PS::Node*>& getNodes() const;

	bool isEventPossible() const;

	std::string toString() const;

private:

	eventType_t eventType;
	double time;
	std::vector<PS::Node*> eventNodes;
	bool isImpossible;

};

} /* namespace Scheduler */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SCHEDULER_EVENT_H_ */
