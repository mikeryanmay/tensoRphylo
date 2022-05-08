/*
 * Edge.cpp
 *
 *  Created on: Aug 22, 2019
 *      Author: xaviermeyer
 */

#include "Edge.h"

#include <cassert>
#include <cmath>
#include <sstream>
#include "Node.h"

namespace Phylogeny {
namespace Structure {


Edge::Edge(Node* aParent, Node *aChild) :
		 parent(aParent), child(aChild) {

	assert(parent != NULL);
	assert(child != NULL);

	// Negative branchs length is node a thing
	// but it might happen during the tree creation from newick string
	length = fabs(parent->age-child->age);
}

Edge::~Edge() {
}


size_t Edge::getId() const {
	assert(child != NULL);
	return child->getId();
}

Node* Edge::getParent() {
	return parent;
}

Node* Edge::getChild(){
	return child;
}


double Edge::getLength() const {
	return length;
}

void Edge::scaleLength(double aFactor) {
	length *=  aFactor;
}

std::string Edge::toString() const {
	std::stringstream ss;
	ss << "Edge [ " << getId() << " ] ( " << length << " ) : " << parent->getId() << " -> " << child->getId() << std::endl;
	return ss.str();
}



} /* namespace Structure */
} /* namespace Phylogeny */
