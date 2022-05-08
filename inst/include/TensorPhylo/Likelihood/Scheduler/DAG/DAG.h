/*
 * DAG.h
 *
 *  Created on: Apr 13, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_DAG_H_
#define LIKELIHOOD_SCHEDULER_DAG_H_


#include <sstream>
#include <vector>
#include <unordered_map>

#include "../IncFwdScheduler.h"

namespace Likelihood {
namespace Scheduler {
namespace DAG {

class NodeDAG;

class DAG {
public:
	DAG(Scheduler::SchedulerSharedPtr aPtrScheduler);
	~DAG();

	void build();
	void rebuild();
	void clear();

	std::string toString() const;

	NodeDAG* getRoot();
	const std::vector<NodeDAG*>& getTips() const;

private:

	Scheduler::SchedulerSharedPtr ptrScheduler;

	NodeDAG* root;
	std::vector<NodeDAG*> tips;

	typedef std::unordered_map<size_t, NodeDAG*> mappingEdgeToDAG_t;
	typedef mappingEdgeToDAG_t::iterator itMapping_t;
	typedef mappingEdgeToDAG_t::const_iterator cItMapping_t;


	void initTips(Event* currentEvent, mappingEdgeToDAG_t &edgeToDAG);
	void asynchronousEvents(Event* currentEvent, mappingEdgeToDAG_t &edgeToDAG);
	void synchronousEvents(Event* currentEvent, mappingEdgeToDAG_t &edgeToDAG);

	NodeDAG* addSegmentedIntegrationNodes(double startTime, double endTime, size_t edgeId, NodeDAG *prevNode);

	void recursiveToString(const NodeDAG* aNode, std::stringstream &ss) const;

};

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SCHEDULER_DAG_H_ */
