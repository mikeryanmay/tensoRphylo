/*
 * Tree.h
 *
 *  Created on: Aug 22, 2019
 *      Author: xaviermeyer
 */

#ifndef DATA_STRUCTURE_TREE_H_
#define DATA_STRUCTURE_TREE_H_

#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <map>

#include "Data/Reader/IncFwdPhyloReader.h"

namespace Phylogeny {
namespace Structure {

class Edge;
class Node;

class Tree {
public:

public:
	Tree(Phylogeny::NexusReader::NexusParserSharedPtr aPtrNexusParser);
	Tree(Phylogeny::NewickReader::NewickParserSharedPtr  aPtrNewickParser);
	~Tree();

	const std::vector<Node*>& getTerminalNodes() const;
	const std::vector<Node*>& getExtinctNodes() const;
	const std::vector<Node*>& getExtantNodes() const;
	const std::vector<Node*>& getNodes() const;

	const std::vector<Edge*>& getEdges() const;

	const Node* getOldestNode() const;

	std::string toString() const;

	std::string getNewickString() const;
	std::string getNewickString(const std::map< int, std::string > &taxaNames) const;

	double getLogLikelihoodCorrection() const;

	void scaleTree(double aFactor);

private:

	static const double LOGOF2;
	double logLikelihoodCorrection;
	Node* oldestNode;

	// Keeping track of different type of nodes
	std::vector<Node*> terminalNodes;
	std::vector<Node*> sampledAncestorNodes;
	std::vector<Node*> extinctNodes;
	std::vector<Node*> extantNodes;
	// All nodes are in nodes
	std::vector<Node*> nodes;

	std::vector<Edge*> edges; // branches

	void buildSingleBranchTree(Phylogeny::NexusReader::NexusParserSharedPtr aPtrNexusParser);
	void buildTreeFromNexus(Phylogeny::NexusReader::NexusParserSharedPtr aPtrNexusParser);
	double createRecursiveNexus(const NxsSimpleNode *nexusNode, Node* treeNode);
	void recursiveScaleTree(double aFactor, Node* treeNode);

	void buildTreeFromNewick(const Phylogeny::NewickReader::TreeNode* aNewickRoot);
	double createRecursiveNewick(const NewickReader::TreeNode *newickNode, Node* treeNode);

	void computeLogLikCorrection();

};

} /* namespace Structure */
} /* namespace Phylogeny */

#endif /* DATA_STRUCTURE_TREE_H_ */
