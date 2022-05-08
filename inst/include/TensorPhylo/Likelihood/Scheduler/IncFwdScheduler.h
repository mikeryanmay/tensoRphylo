/*
 * IncFwdScheduler.h
 *
 *  Created on: Aug 29, 2019
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_INCFWDSCHEDULER_H_
#define LIKELIHOOD_SCHEDULER_INCFWDSCHEDULER_H_

#include <boost/smart_ptr/shared_ptr.hpp>
#include "EventType.h"

namespace Likelihood {

namespace Scheduler {
	class BaseScheduler;
	class Event;
	typedef boost::shared_ptr<BaseScheduler> SchedulerSharedPtr;
}
}

#endif /* LIKELIHOOD_SCHEDULER_INCFWDSCHEDULER_H_ */
