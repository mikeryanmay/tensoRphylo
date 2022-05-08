//    HOGAN is an implementation of a parallel Metropolis-Hastings algorithm 
//    developped for evolutionnary biology model.
//    Copyright (C) 2016  Xavier Meyer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
/**
 * @file TreeNode.h
 *
 * @date Jan 21, 2015
 * @author meyerx
 * @brief
 */
#ifndef TREENODEDL_H_
#define TREENODEDL_H_

#include <boost/fusion/adapted/adt/adapt_adt.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <stddef.h>
#include <sstream>
#include <vector>
#include <map>

#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <boost/fusion/adapted/struct/detail/adapt_base_attr_filler.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/comparison/not_equal.hpp>
#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/logical/bool.hpp>
#include <boost/preprocessor/repetition/detail/for.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/tuple/eat.hpp>
#include <boost/preprocessor/tuple/elem.hpp>

namespace Phylogeny {
namespace NewickReader {

class TreeNode;

typedef std::vector<TreeNode> TreeNode_children;

class TreeNode {
public:

	TreeNode();
	~TreeNode();

	size_t getId() const;
	std::string getName() const;
	int getIdTaxa() const;
	double getLength() const;
	const TreeNode_children& getChildren() const;

	bool hasMetadata() const;
	std::map<std::string, std::string> extractMetadata() const;
	size_t defineRBCompatibleNodeID() const;

	bool isLeaf() const;

	std::string toString() const;
	std::string subtreeToString() const;

	static void resetIdSeq();

private:
	static size_t idSeq;
	size_t id;

public:
	int idTaxa;
	std::string name;
	std::vector<std::string> metadata;
	double length;
	TreeNode_children children;

	void memorizeBLRecursively(std::map<size_t, double> &branchLengths) const;
	void memorizeNamesRecursively(std::vector<std::string> &taxaNames);

};

} /* namespace NewickReader */
} /* namespace Phylogeny */

BOOST_FUSION_ADAPT_STRUCT(Phylogeny::NewickReader::TreeNode,
		(Phylogeny::NewickReader::TreeNode_children, children)
		(std::string, name)
		(std::vector<std::string>, metadata)
		(double, length) )


#endif /* TREENODEDL_H_ */
