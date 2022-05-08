/*
 * Node.h
 *
 *  Created on: Aug 22, 2019
 *      Author: xaviermeyer
 */

#ifndef DATA_STRUCTURE_NODE_H_
#define DATA_STRUCTURE_NODE_H_

#include <map>
#include <string>
#include <vector>

namespace Phylogeny {
namespace Structure {

class Edge;

class Node {
public:
	Node(double aAge, int aIdTaxa);
	Node(size_t aId, double aAge,  int aIdTaxa);
	~Node();

	const std::vector<Edge*>&  getEdges() const;

	Edge* getEdgeToParent() const;
	std::vector<Edge*> getEdgesToChildren() const;

	Node* getParentNode() const;
	std::vector<Node*> getChildrenNodes() const;

	void defineAsRootNode(bool aIsRoot);
	void defineAsOriginNode(bool aIsOrigin);
	void defineAsGhostNode(bool aIsGhost);

	bool isOriginNode() const;
	bool isRootNode() const;
	bool isTerminalNode() const;
	bool isGhostNode() const;

	bool isExtinct() const; // Is terminal and has age>0
	bool isExtant() const; // Is terminal and has age=0
	bool isSpeciationNode() const; // Not Extinct or extant and as two descendant and one ancestor
	bool isSampledAncestor() const; // Is an internal node with one descendant

	size_t getId() const;
	int getTaxaId() const;
	std::string getName(const std::vector<std::string> &aTaxaNames) const;

	void setAge(double aAge);
	double getAge() const;

	void addEdge(Edge* aEdge);

	std::string toString() const;
	std::string buildSubtreeNewickString(const std::map<int, std::string> &taxaIdToNameMap) const;

	double getSubtreeBranchLength() const;

protected:
	static size_t idSeq;

	size_t id;
	int idTaxa;
	double age;
	bool isOrigin, isRoot, isGhost;

	std::vector<Edge*> edges;

private:

	friend class Edge;

	void addNodeToOrderedString(std::string &str, const Node *parent, const std::map<int, std::string> &taxaIdToNameMap) const;

	double recursiveBranchLengthSum(const Node* aNode) const;

};

} /* namespace Structure */
} /* namespace Phylogeny */

#endif /* DATA_STRUCTURE_NODE_H_ */
