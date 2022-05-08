
/**
 * @file TreeNode.cpp
 *
 * @date Jan 22, 2015
 * @author meyerx
 * @brief
 */

#include "TreeNode.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cassert>
#include <iostream>

namespace Phylogeny {
namespace NewickReader {

size_t TreeNode::idSeq = 0;

TreeNode::TreeNode() : id(idSeq++), idTaxa(-1), length(0.) {

}


TreeNode::~TreeNode() {

}

std::string TreeNode::toString() const {
	using std::stringstream;

	stringstream ss;
	ss << "Node [" << id << "]" << std::endl << "Name : " << name << std::endl << "Id Taxa : " << idTaxa << std::endl << "Length : " << length << std::endl;
	//assert(metadata.size() == values.size());
	for(size_t i=0; i<metadata.size(); ++i) {
		ss << "Metadata : " << metadata[i] << std::endl;// " : " << values[i] << std::endl;
	}

	for(size_t i=0; i<children.size(); ++i){
		ss << "[" << children[i].id << "]" << children[i].name;
		if(i < children.size()-1) ss << " :: ";
	}

	return ss.str();

}

std::string TreeNode::subtreeToString() const {
	std::stringstream ss;
	ss << "-------------------------------------------" << std::endl;
	ss << toString() << std::endl;

	for (size_t i = 0; i < children.size(); ++i) {
		ss << children[i].subtreeToString();
	}

	return ss.str();
}

void TreeNode::resetIdSeq() {
	idSeq = 0;
}

size_t TreeNode::getId() const {
	return id;
}

std::string TreeNode::getName() const {
	return name;
}

int TreeNode::getIdTaxa() const {
	return idTaxa;
}

double TreeNode::getLength() const {
	return length;
}

const TreeNode_children& TreeNode::getChildren() const {
	return children;
}


bool TreeNode::hasMetadata() const {
	return !metadata.empty();
}

std::map<std::string, std::string> TreeNode::extractMetadata() const {

	std::map<std::string, std::string> metadataMap;
	for(size_t iM=0; iM<metadata.size();++iM) {
	    std::vector<std::string> result;
	    boost::algorithm::split(result, metadata[iM], boost::algorithm::is_any_of("="));
	    assert(result.size() == 2);
	    metadataMap[result[0]] = result[1];
	}

	return metadataMap;
}

bool TreeNode::isLeaf() const {
	return children.empty();
}

/**
 * This function extract the compatible revbayes node id if possible.
 * If not it uses the default newick node id.
 */
size_t TreeNode::defineRBCompatibleNodeID() const {
	size_t idNode = this->getId();
	if(this->hasMetadata()) {
		std::map<std::string, std::string> metadataMap = this->extractMetadata();
		std::map<std::string, std::string>::iterator itFind = metadataMap.find("index");
		if(itFind != metadataMap.end()) {
			idNode = std::stol(itFind->second);
		}
	}
	return idNode;
}

void TreeNode::memorizeBLRecursively(std::map<size_t, double> &branchLengths) const {

	branchLengths.insert(std::make_pair(id, length));

	// For all children do the same
	for(size_t iC=0; iC<children.size(); ++iC) {
		children[iC].memorizeBLRecursively(branchLengths);
	}
}


void TreeNode::memorizeNamesRecursively(std::vector<std::string> &taxaNames) {

	if(name.empty()) {
		idTaxa = -1;
	} else {
		idTaxa = taxaNames.size();
		taxaNames.push_back(name);
	}

	// For all children do the same
	for(size_t iC=0; iC<children.size(); ++iC) {
		children[iC].memorizeNamesRecursively(taxaNames);
	}
}

} /* namespace NewickReader */
} /* namespace Phylogeny */
