/*
 * Edge.h
 *
 *  Created on: Aug 22, 2019
 *      Author: xaviermeyer
 */

#ifndef DATA_STRUCTURE_EDGE_H_
#define DATA_STRUCTURE_EDGE_H_

#include <iostream>

namespace Phylogeny {
namespace Structure {

class Node;

class Edge {
public:
	Edge(Node* aParent, Node *aChild);
	~Edge();

	size_t getId() const;
	Node* getParent();
	Node* getChild();
	double getLength() const;

	void scaleLength(double aFactor);

	std::string toString() const;

protected:

	Node *parent, *child;
	double length;

private:

	friend class Node;
};

} /* namespace Structure */
} /* namespace Phylogeny */

#endif /* DATA_STRUCTURE_EDGE_H_ */
