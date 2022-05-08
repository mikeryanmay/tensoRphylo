/*
 * Tree.cpp
 *
 *  Created on: Aug 22, 2019
 *      Author: xaviermeyer
 */

#include "Tree.h"

#include <iostream>
#include <limits>
#include <sstream>
#include <limits>

#include "Edge.h"
#include "Node.h"

#include "Data/Reader/IncPhyloReader.h"

namespace Phylogeny {
namespace Structure {

const double Tree::LOGOF2 = log(2.);


Tree::Tree(Phylogeny::NexusReader::NexusParserSharedPtr aPtrNexusParser) {

	if(aPtrNexusParser->isSingleBranchHack()) {
		buildSingleBranchTree(aPtrNexusParser);
	} else {
		buildTreeFromNexus(aPtrNexusParser);
	}
	computeLogLikCorrection();
}

Tree::Tree(Phylogeny::NewickReader::NewickParserSharedPtr aPtrNewickParser) {
	buildTreeFromNewick(&aPtrNewickParser->getRoot());
	computeLogLikCorrection();
}

Tree::~Tree() {
	for(size_t iE=0; iE<edges.size(); ++iE) {
		delete edges[iE];
	}

	for(size_t iN=0; iN<nodes.size(); ++iN) {
		delete nodes[iN];
	}
}

void Tree::buildSingleBranchTree(Phylogeny::NexusReader::NexusParserSharedPtr aPtrNexusParser) {

	const std::vector<double> &nodeTimes = aPtrNexusParser->getNodeTimes();
	assert(nodeTimes.size() >= 2 && "There should be at least a child and a parent node to define a single branch.");

	// Create first node
	Node *childNode = new Node(0., 0);
	nodes.push_back(childNode);
	assert(nodeTimes[0] == 0. && "The first node of a single branch should always be at time 0.");

	const std::vector<double> &pVec = aPtrNexusParser->getProbForTaxaId(0);
	bool isGhostNode = true;
	for(size_t iP=0; iP<pVec.size(); ++iP) {
		isGhostNode = isGhostNode && pVec[iP] == 0.;
	}
	childNode->defineAsGhostNode(isGhostNode);

	extantNodes.push_back(childNode);
	terminalNodes.push_back(childNode);

	for(size_t iN=1; iN < nodeTimes.size(); ++iN) {

		Node *parentNode;
		if(iN < nodeTimes.size()-1) parentNode = new Node(nodeTimes[iN], iN);
		else parentNode = new Node(nodeTimes[iN], -1);

		nodes.push_back(parentNode);
		Edge* edge = new Edge(parentNode, childNode);
		edges.push_back(edge);

		childNode->addEdge(edge);
		parentNode->addEdge(edge);

		childNode = parentNode;
	}

	// Origination node is now in childNode
	childNode->defineAsRootNode(true);
	childNode->defineAsOriginNode(true);

	// set the oldest node
	oldestNode = childNode;
}

void Tree::buildTreeFromNexus(Phylogeny::NexusReader::NexusParserSharedPtr aPtrNexusParser) {

	const NxsSimpleTree *nxsTree = aPtrNexusParser->getNxsTree();

	int idTaxa = nxsTree->GetRootConst()->GetTaxonIndex() == std::numeric_limits<unsigned int>::max() ? -1 : nxsTree->GetRootConst()->GetTaxonIndex();
	Node *root = new Node(0., idTaxa);
	nodes.push_back(root);	// Register root node
	root->defineAsRootNode(true);

	// Deal with origination node
	double originTimeRespToRoot = nxsTree->GetRootConst()->GetEdgeToParent().GetDblEdgeLen();
	if(originTimeRespToRoot > 0.) {
		oldestNode = new Node(0., -1); 							// Origination
		nodes.push_back(oldestNode);							// Register origination node
		root->setAge(originTimeRespToRoot);					// Set root age respective to origination
		Edge* edge = new Edge(oldestNode, root);		// Create edge
		edges.push_back(edge);										// Register edge
		oldestNode->addEdge(edge);								// Add edge to child
		root->addEdge(edge);											// Add edge to parent
	} else {
		oldestNode = root;
	}
	oldestNode->defineAsOriginNode(true);

	double maxTempAge = createRecursiveNexus(nxsTree->GetRootConst(), root);

	// The current node ages are from past (0) to present (maxTempAge)
	// Branch length have the absolute age different, so they are correct.
	// We need however to re-date each node
	for(size_t iN=0; iN<nodes.size(); ++iN) {
		double oldAge = nodes[iN]->getAge();
		double newAge = maxTempAge - oldAge;

		if(nodes[iN]->isTerminalNode() && newAge < 1.e-4) {
			newAge = 0.;
		}

		nodes[iN]->setAge(newAge);
	}

	// Keep track of sampled ancestor, fossil and leaf nodes
	for(size_t iN=0; iN<nodes.size(); ++iN) {
		if(nodes[iN]->isSampledAncestor()) {
			sampledAncestorNodes.push_back(nodes[iN]);
		} else if(nodes[iN]->isExtant()) {
			extantNodes.push_back(nodes[iN]);
		} else if(nodes[iN]->isExtinct()) {
			extinctNodes.push_back(nodes[iN]);
		}
	}
}

double Tree::createRecursiveNexus(const NxsSimpleNode *nexusNode, Node* treeNode) {

	double maxAge = treeNode->getAge();

	// For each child of newick Node, we create the new treeNode
	for(size_t iC=0; iC<nexusNode->GetChildren().size(); ++iC) {
		const NxsSimpleNode *childNexus= nexusNode->GetChildren()[iC];
		double branchLengthToParent = childNexus->GetEdgeToParent().GetDblEdgeLen();

		// Creating new node and keeping track of oldest age (with root having 0)
		double temporaryAge = treeNode->getAge() + branchLengthToParent;


		int idTaxa = childNexus->GetTaxonIndex() == std::numeric_limits<unsigned int>::max() ? -1 : childNexus->GetTaxonIndex();
		Node* childNode = new Node(temporaryAge, idTaxa);
		double tempMaxAge = createRecursiveNexus(childNexus, childNode);
		maxAge = std::max(tempMaxAge, maxAge);

		// Updating tree structure
		nodes.push_back(childNode);									// Register node
		Edge* edge = new Edge(treeNode, childNode);		// Create edge
		edges.push_back(edge);											// Register edge
		childNode->addEdge(edge);										// Add edge to child
		treeNode->addEdge(edge);										// Add edge to parent
	}

	assert(treeNode);
	if(treeNode->isTerminalNode()) {									// Register terminal nodes
		terminalNodes.push_back(treeNode);
	}

	return maxAge;

}

void Tree::buildTreeFromNewick(const NewickReader::TreeNode* aNewickRoot) {

	size_t nodeId = aNewickRoot->defineRBCompatibleNodeID();
	Node *root = new Node(nodeId, 0., aNewickRoot->getIdTaxa());
	nodes.push_back(root);	// Register root node
	root->defineAsRootNode(true);

	// Deal with origination node
	double originTimeRespToRoot = aNewickRoot->getLength();
	if(originTimeRespToRoot > 0.) {
		oldestNode = new Node(0., -1); 							// Origination
		nodes.push_back(oldestNode);							// Register origination node
		root->setAge(originTimeRespToRoot);					// Set root age respective to origination
		Edge* edge = new Edge(oldestNode, root);		// Create edge
		edges.push_back(edge);										// Register edge
		oldestNode->addEdge(edge);								// Add edge to child
		root->addEdge(edge);											// Add edge to parent
	} else {
		oldestNode = root;
	}
	oldestNode->defineAsOriginNode(true);

	double maxTempAge = createRecursiveNewick(aNewickRoot, root);

	// The current node ages are from past (0) to present (maxTempAge)
	// Branch length have the absolute age different, so they are correct.
	// We need however to re-date each node
	for(size_t iN=0; iN<nodes.size(); ++iN) {
		double oldAge = nodes[iN]->getAge();
		double newAge = maxTempAge - oldAge;

		if(nodes[iN]->isTerminalNode() && newAge < 1.e-4) {
			newAge = 0.;
		}

		nodes[iN]->setAge(newAge);
	}

	// Keep track of fossil and leaf nodes
	for(size_t iN=0; iN<terminalNodes.size(); ++iN) {
		if(terminalNodes[iN]->isExtant()) {
			extantNodes.push_back(terminalNodes[iN]);
		} else {
			assert(terminalNodes[iN]->isExtinct());
			extinctNodes.push_back(terminalNodes[iN]);
		}
	}
}

double Tree::createRecursiveNewick(const NewickReader::TreeNode *newickNode, Node* treeNode) {

	double maxAge = treeNode->getAge();

	// For each child of newick Node, we create the new treeNode
	for(size_t iC=0; iC<newickNode->getChildren().size(); ++iC) {
		const NewickReader::TreeNode *childNewick =&newickNode->getChildren()[iC];
		double branchLengthToParent = childNewick->getLength();

		// Creating new node and keeping track of oldest age (with root having 0)
		double temporaryAge = treeNode->getAge() + branchLengthToParent;

		size_t nodeId = childNewick->defineRBCompatibleNodeID();
		Node* childNode = new Node(nodeId, temporaryAge, childNewick->getIdTaxa());
		double tempMaxAge =  createRecursiveNewick(childNewick, childNode);
		maxAge = std::max(tempMaxAge, maxAge);

		// Updating tree structure
		nodes.push_back(childNode);									// Register node
		Edge* edge = new Edge(treeNode, childNode);		// Create edge
		edges.push_back(edge);											// Register edge
		childNode->addEdge(edge);										// Add edge to child
		treeNode->addEdge(edge);										// Add edge to parent
	}

	assert(treeNode);
	if(treeNode->isTerminalNode()) {									// Register terminal nodes
		terminalNodes.push_back(treeNode);
	}

	return maxAge;

}

const std::vector<Node*>& Tree::getTerminalNodes() const {
	return terminalNodes;
}

const std::vector<Node*>& Tree::getExtantNodes() const {
	return extantNodes;
}

const std::vector<Node*>& Tree::getExtinctNodes() const {
	return extinctNodes;
}

const std::vector<Node*>& Tree::getNodes() const {
	return nodes;
}

const std::vector<Edge*>& Tree::getEdges() const {
	return edges;
}

const Node* Tree::getOldestNode() const {
	return oldestNode;
}

std::string Tree::toString() const {
	std::stringstream ss;
	for(size_t iN=0; iN<nodes.size(); ++iN) {
		ss << nodes[iN]->toString() << std::endl;
	}
	return ss.str();
}


std::string Tree::getNewickString() const {
	std::map<int, std::string> taxaIdToNameMap;
	return oldestNode->buildSubtreeNewickString(taxaIdToNameMap);
}

std::string Tree::getNewickString(const std::map<int, std::string> &taxaIdToNameMap) const {
	return oldestNode->buildSubtreeNewickString(taxaIdToNameMap);
}

double Tree::getLogLikelihoodCorrection() const {
	return logLikelihoodCorrection;
}

void Tree::computeLogLikCorrection(){

	size_t nSampledAncestor = sampledAncestorNodes.size();
	size_t nExtinct = extinctNodes.size();
	size_t nExtant = extantNodes.size();
	size_t nTaxa = nSampledAncestor + nExtinct + nExtant;

	double logFactorial = 0.;
	for(size_t iF=0; iF<(nTaxa-nExtinct); ++iF) {
		logFactorial += log(iF+1);
	}

	logLikelihoodCorrection = LOGOF2*(nTaxa - nSampledAncestor - 1.0);
	logLikelihoodCorrection -= logFactorial;
}


void Tree::scaleTree(double aFactor) {
	recursiveScaleTree(aFactor, oldestNode);
}

void Tree::recursiveScaleTree(double aFactor, Node* treeNode) {

		double newNodeAge = aFactor*treeNode->getAge();
		treeNode->setAge(newNodeAge);

		std::vector<Edge*> edgesToChildren = treeNode->getEdgesToChildren();
		for(size_t iE=0; iE<edgesToChildren.size(); ++iE) {
			edgesToChildren[iE]->scaleLength(aFactor);
			recursiveScaleTree(aFactor, edgesToChildren[iE]->getChild());
		}
}

} /* namespace Structure */
} /* namespace Phylogeny */
