/*
 * Node.cpp
 *
 *  Created on: Aug 22, 2019
 *      Author: xaviermeyer
 */

#include "Node.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

#include "Edge.h"

namespace Phylogeny {
namespace Structure {

size_t Node::idSeq = 0;

Node::Node(double aAge, int aIdTaxa) : id(idSeq++), idTaxa(aIdTaxa), age(aAge) {
	isOrigin = false;
	isRoot = false;
	isGhost = false;
}

Node::Node(size_t aId, double aAge,  int aIdTaxa) : id(aId), idTaxa(aIdTaxa), age(aAge){
	isOrigin = false;
	isRoot = false;
	isGhost = false;
}

Node::~Node() {
}


const std::vector<Edge*>&  Node::getEdges() const {
	return edges;
}

Edge*  Node::getEdgeToParent() const {
	Edge* edgeToParent = NULL;

	for(size_t iE=0; iE<edges.size(); ++iE) {
		if(edges[iE]->child == this) {
			edgeToParent = edges[iE];
		}
	}
	return edgeToParent;
}

std::vector<Edge*> Node::getEdgesToChildren() const {
	std::vector<Edge*> edgesToChildren;

	for(size_t iE=0; iE<edges.size(); ++iE) {
		if(edges[iE]->parent == this) {
			edgesToChildren.push_back(edges[iE]);
		}
	}
	return edgesToChildren;
}

Node* Node::getParentNode() const {
	Edge* edgeToParent  = this->getEdgeToParent();
	return edgeToParent->parent;
}

std::vector<Node*> Node::getChildrenNodes() const {
	std::vector<Edge*> edgesToChildren = this->getEdgesToChildren();

	std::vector<Node*> childrenNodes;
	for(size_t iE = 0; iE < edgesToChildren.size(); ++iE) {
		childrenNodes.push_back(edgesToChildren[iE]->child);
	}

	return childrenNodes;
}

void Node::defineAsRootNode(bool aIsRoot) {
	isRoot = aIsRoot;
}

void Node::defineAsOriginNode(bool aIsOrigin) {
	isOrigin = aIsOrigin;
}

void Node::defineAsGhostNode(bool aIsGhost) {
	isGhost = aIsGhost;
}

bool Node::isOriginNode() const {
	return isOrigin;
}

bool Node::isRootNode() const {
	return isRoot;
}

bool Node::isTerminalNode() const {
	// One could look if edges.size() == 1, but less safe
	std::vector<Node*> childrenNodes;
	childrenNodes = this->getChildrenNodes();
	return childrenNodes.empty();
}

bool Node::isGhostNode() const {
	return isGhost;
}


bool Node::isExtinct() const {
	return this->isTerminalNode() && age > 0;
}

bool Node::isSpeciationNode() const {

	return !this->isExtant() && !this->isExtinct() && age > 0 &&
			this->getEdgesToChildren().size() == 2;
}

bool Node::isExtant() const {
	return this->isTerminalNode() && age == 0.;
}

bool Node::isSampledAncestor() const {
	return this->getEdgesToChildren().size() == 1 && isOrigin == false;
}


size_t Node::getId() const {
	return id;
}


int Node::getTaxaId() const {
	return idTaxa;
}

std::string Node::getName(const std::vector<std::string> &aTaxaNames) const {
	if(idTaxa >= 0) {
		assert(idTaxa < (int)aTaxaNames.size());
		return aTaxaNames[idTaxa];
	} else {
		return std::string();
	}
}


void Node::setAge(double aAge) {
	if(isTerminalNode() && fabs(aAge) < 2.*std::numeric_limits<double>::epsilon()) {
		age = 0;
	} else {
		assert(aAge >= 0.);
		age = aAge;
	}
}

double Node::getAge() const {
	return age;
}

void Node::addEdge(Edge* aEdge) {
	assert(aEdge->parent == this || aEdge->child == this); // Check that this is a meaningful edges
	for(size_t iE=0; iE<edges.size(); ++iE) { // Check that each edge is unique
		assert(!(aEdge->child == edges[iE]->child && aEdge->parent == edges[iE]->parent));
	}
	edges.push_back(aEdge);
}

std::string Node::toString() const {
	std::stringstream ss;
	ss << "Node [" << id <<  " : " << idTaxa << "] : " << " - " << std::scientific << age << std::endl;
	for(size_t iE=0; iE<edges.size(); ++iE) {
		ss << " -- Edge : " << edges[iE]->getParent()->getId() << " --> " << edges[iE]->getChild()->getId() << std::endl;
	}

	/*std::cout << (isExtant() ? "extant" : "not extant" ) << std::endl;
	std::cout << (isExtinct() ? "extinct" : "not extinct" ) << std::endl;
	std::cout << (isTerminalNode() ? "Terminal" : "not terminal" ) << std::endl;
	std::cout << (age == 0. ? "present" : "not present" ) << std::endl;
	std::cout << (isSpeciationNode() ? "speciation" : "not speciation" ) << std::endl;
	std::cout << (isSampledAncestor() ? "sampledAnc" : "not sampledAnc" ) << std::endl;*/

	return ss.str();
}

std::string Node::buildSubtreeNewickString(const std::map<int, std::string> &taxaIdToNameMap) const {
	std::vector<double> branchLength;
	std::string nwkStr;
	addNodeToOrderedString(nwkStr, NULL, taxaIdToNameMap);
	nwkStr.append(";");
	return nwkStr;
}

double Node::getSubtreeBranchLength() const {
	return recursiveBranchLengthSum(this);
}

void Node::addNodeToOrderedString(std::string &str, const Node *parent,
		const std::map<int, std::string> &taxaIdToNameMap) const {

	if(!this->isTerminalNode()) {
		std::vector<Node*> children = this->getChildrenNodes();

		str.push_back('(');
		for(size_t iC=0; iC<children.size(); ++iC) {
			assert(children[iC]);
			children[iC]->addNodeToOrderedString(str, this, taxaIdToNameMap);
			if(iC<children.size()-1) str.push_back(',');
		}
		str.push_back(')');
	} else {
		assert(idTaxa >= 0);
	}

	std::stringstream ss;
	if(idTaxa >= 0) {
		if(taxaIdToNameMap.empty()) {
			ss << idTaxa;
		} else {
			std::map<int, std::string>::const_iterator it = taxaIdToNameMap.find(idTaxa);
			assert(it != taxaIdToNameMap.end());
			ss << it->second;
		}
	}
	if(this->getEdgeToParent() != NULL) {
		ss << ":" << this->getEdgeToParent()->getLength();
	}
	str.append(ss.str());

}

double Node::recursiveBranchLengthSum(const Node* aNode) const {
	double sumBL = 0.;

	std::vector<Edge*> edgesToChildren = getEdgesToChildren();
	for(size_t i=0; i<edgesToChildren.size(); ++i) {
		sumBL += recursiveBranchLengthSum(edgesToChildren[i]->getChild());
		sumBL += edgesToChildren[i]->getLength();
	}

	return sumBL;
}

} /* namespace Structure */
} /* namespace Phylogeny */
