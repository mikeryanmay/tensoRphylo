/*
 * SchedulerDAG.h
 *
 *  Created on: Apr 14, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_DAG_SCHEDULERDAG_H_
#define LIKELIHOOD_SCHEDULER_DAG_SCHEDULERDAG_H_

#include <list>
#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>

namespace Likelihood {
namespace Scheduler {
namespace DAG {

class NodeDAG;

class DAG;
typedef boost::shared_ptr<DAG> DAGSharedPtr;

class SchedulerDAG {
public:
	SchedulerDAG(DAGSharedPtr aPtrDAG);
	~SchedulerDAG();

	bool isDone();
	void reset();

	NodeDAG* getNextTask();

	void signalDone(NodeDAG* aNode);

private:

	DAGSharedPtr ptrDAG;

	bool rootSignaled;
	std::list<NodeDAG*> pendingTasks;

	void init();

};

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SCHEDULER_DAG_SCHEDULERDAG_H_ */
