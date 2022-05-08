/*
 * NodeDAG.h
 *
 *  Created on: Apr 13, 2020
 *      Author: xaviermeyer
 */

#ifndef LIKELIHOOD_SCHEDULER_DAG_NODEDAG_H_
#define LIKELIHOOD_SCHEDULER_DAG_NODEDAG_H_

#include "NodeDAG.h"
#include <set>
#include <string>
#include <vector>

namespace Likelihood {

namespace StateType {
namespace Vector {

class EigenState; // Fwd declaration
// TODO replace by a template version of the DAG suite for any type of State.

}
}

namespace Scheduler {

class Event;

namespace DAG {

typedef Likelihood::StateType::Vector::EigenState StateType;

typedef enum {
	INTEGRATE_TASK,
	EVENT_TASK
} taskType_t;

class NodeDAG {
public:

	NodeDAG(double aTime, Event *event, long int aIdEdge);
	NodeDAG(double aStartTime, double aEndType, long int aIdEdge);
	~NodeDAG();

	void attachState(StateType *aPState);
	void detachState(StateType *aPState);
	void detachStates();

	std::vector<StateType*>& getProbStates();
	const std::vector<StateType*>& getProbStates() const;

	bool isReady() const;
	bool isIntegrationTask() const;
	bool isEventTask() const;

	size_t getId() const;
	double getStartTime() const;
	double getEndTime() const;
	double getTime() const;

	Event* getEvent();

	void addChildNode(NodeDAG* aNode);
	void addParentNode(NodeDAG* aNode);

	void removeChildNode(NodeDAG* aNode);
	void removeParentNode(NodeDAG* aNode);

	long int getEdgeIdMapping() const;

	void setChildNodeReady(NodeDAG* aChildNode);
	void signalReadyToParents();
	void resetChildReadyState();

	std::vector<NodeDAG*>& getChildNodes();
	const std::vector<NodeDAG*>& getChildNodes() const;

	std::vector<NodeDAG*>& getParentNodes();
	const std::vector<NodeDAG*>& getParentNodes() const;

	std::string toString() const;

	void identifyTips(std::vector<NodeDAG*> &tips, std::set<size_t> &visitedTips);
	void recursiveClear();

	double identifySubDAGSumBL();

private:

	static size_t idSeq;

	size_t id;
	taskType_t taskType;
	Event *event;

	long int idIntEdge;
	double startTime, endTime;

	typedef std::vector<StateType*> vecStates_t;
	typedef vecStates_t::iterator itVecStates_t;
	typedef vecStates_t::const_iterator cItVecStates_t;

	vecStates_t probStates;

	std::vector<NodeDAG*> childNodes;
	std::vector<bool> childNodesReady;
	std::vector<NodeDAG*> parentNodes;


	double identifySubDAGSumBL(NodeDAG *aNode);

};

} /* namespace DAG */
} /* namespace Scheduler */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_SCHEDULER_DAG_NODEDAG_H_ */
