/*
 * ParallelSchedulerDAG.h
 *
 *  Created on: Apr 14, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_DAG_PARALLELSCHEDULERDAG_H_
#define LIKELIHOOD_SCHEDULER_DAG_PARALLELSCHEDULERDAG_H_


#include "Utils/Parallel/Manager.h"
#if defined(_OPENMP)

#include <list>
#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <omp.h>

namespace Likelihood {
namespace Scheduler {
namespace DAG {

class NodeDAG;

class DAG;
typedef boost::shared_ptr<DAG> DAGSharedPtr;

class ParallelSchedulerDAG {
public:
	ParallelSchedulerDAG(size_t nThread, DAGSharedPtr aPtrDAG);
	~ParallelSchedulerDAG();

	bool isDone();
	void reset();

	NodeDAG* getNextTask(size_t iThread);

	bool tryNonBlockingSignalAndPop(size_t iThread, NodeDAG* &inOutNode);

	void signalDone(size_t iThread, NodeDAG* aNode);

	size_t getNThread() const;

private:

	const size_t N_THREADS;
	DAGSharedPtr ptrDAG;

	bool rootSignaled;
	std::vector<omp_lock_t> locks;
	std::vector< std::list<NodeDAG*> > pendingTasks;

	void init();

	NodeDAG* getNextTaskFromLocalThreadPool(size_t iThread);
	NodeDAG* getNextTaskFromGlobalThreadPool(size_t iThread);

};

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SCHEDULER_DAG_PARALLELSCHEDULERDAG_H_ */

#endif
