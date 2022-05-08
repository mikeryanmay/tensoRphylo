/*
 * EventType.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_EVENTTYPE_H_
#define LIKELIHOOD_SCHEDULER_EVENTTYPE_H_

namespace Likelihood {
namespace Scheduler {

typedef enum {
	PRESENT_TIME_EVENT,
	NODE_EVENT,
	SYCHRONOUS_SPECIATION_EVENT,
	SYNCHRONOUS_EXTINCTION_EVENT,
	SYNCHRONOUS_SAMPLING_EVENT,
	SYNCHRONOUS_DESTRUCTIVE_SAMPLING_EVENT,
	SYNCHRONOUS_RATE_SHIFT,
	SYNCHRONOUS_RESCALING_EVENT,
	FINAL_NODE_EVENT,
	SYNCHRONOUS_MONITORING_PROB
} eventType_t;

}
}

#endif /* LIKELIHOOD_SCHEDULER_EVENTTYPE_H_ */