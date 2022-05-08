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
 * @file NewickParser.h
 *
 * @date Jan 21, 2015
 * @author meyerx
 * @brief
 */
#ifndef NEWICKPARSER_H_
#define NEWICKPARSER_H_

#include <fstream>
#include <iostream>
#include <string>

#include "NewickGrammar.h"
#include "TreeNode.h"

namespace Phylogeny {
namespace NewickReader {

class NewickParser {
public:
	enum stringType_t {IS_FILE_NAME=0, IS_TREE_STRING=1};
public:
	NewickParser(const std::string &aString, const stringType_t stringType=IS_FILE_NAME);
	~NewickParser();

	const TreeNode& getRoot() const;
	const std::map<size_t, double>& getInitialBLs() const;
	const std::vector<std::string>& getTaxaNames() const;

private:

	const std::string fileName;
	TreeNode root;
	std::vector<std::string> taxaNames;

	std::map<size_t, double> initialBLs;

	bool readFile(std::string &str);
	bool parseFile(std::string &str);

	void defineInitialBranchLengthRecursively();
	void defineNamesRecursively();

};

} /* namespace NewickReader */
} /* namespace Phylogeny */

#endif /* NEWICKPARSER_H_ */
