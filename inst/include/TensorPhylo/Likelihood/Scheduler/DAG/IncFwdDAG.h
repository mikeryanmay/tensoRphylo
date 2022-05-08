/*
 * IncFwdScheduler.h
 *
 *  Created on: Apr 14, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_DAG_INCFWDDAG_H_
#define LIKELIHOOD_SCHEDULER_DAG_INCFWDDAG_H_

#include <boost/smart_ptr/shared_ptr.hpp>

namespace Likelihood {
namespace Scheduler {
namespace DAG {

	class NodeDAG;

	class DAG;
	typedef boost::shared_ptr<DAG> DAGSharedPtr;

	class SchedulerDAG;
	typedef boost::shared_ptr<SchedulerDAG> SchedulerDAGSharedPtr;

	class ParallelSchedulerDAG;
	typedef boost::shared_ptr<ParallelSchedulerDAG> ParallelSchedulerDAGSharedPtr;
}
}
}

#endif /* LIKELIHOOD_SCHEDULER_DAG_INCFWDDAG_H_ */
