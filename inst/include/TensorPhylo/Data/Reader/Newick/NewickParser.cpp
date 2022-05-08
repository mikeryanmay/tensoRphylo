/**
 * @file NewickParser.cpp
 *
 * @date Jan 22, 2019
 * @author meyerx
 * @brief
 */
#include "NewickParser.h"

#include <boost/spirit/home/qi/numeric/detail/real_impl.hpp>
#include <fstream>
#include <iosfwd>

namespace Phylogeny {
namespace NewickReader {

NewickParser::NewickParser(const std::string &aString, const stringType_t stringType) :
		fileName(aString) {

	std::string str;
	if(stringType == IS_FILE_NAME) {
		readFile(str);
	} else {
		str = aString;
	}
	parseFile(str);
	defineNamesRecursively();
	defineInitialBranchLengthRecursively();
}

NewickParser::~NewickParser() {
}

bool NewickParser::readFile(std::string &str) {
	using std::ifstream;

	ifstream iFile(fileName.c_str());
	if (!iFile.good()) {
		std::cerr << "bool NewickParser::readFile(string &str)" << std::endl;
		std::cerr << "Newick file not found." << std::endl;
		return false;
	}

	getline(iFile, str);

	return !str.empty();
}

bool NewickParser::parseFile(std::string &str) {
	namespace qi = ::boost::spirit::qi;

	//str = "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);";
	//str = "(((One:0.1,Two:0.2)Sub1:0.3,(Three:0.4,Four:0.5)Sub2:0.6)Sub3:0.7,Five:0.8)Root:0.9;";
	//str = "(((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;";
	//std::cout << str << std::endl;

	NewickGrammar<std::string::const_iterator> grammar;
	// Parse
	std::string::const_iterator iter = str.begin();
	std::string::const_iterator end = str.end();
	bool result = qi::phrase_parse(iter, end, grammar, qi::space, root);

	//std::cout << root << std::endl;

	if (!result) {
		std::cerr << "bool NewickParser::parseFile(string &str)" << std::endl;
		std::cerr << "Newick file badly formated." << std::endl;
		return false;
	}

	TreeNode::resetIdSeq();

	return true;
}

const TreeNode& NewickParser::getRoot() const {
	return root;
}


const std::vector<std::string>& NewickParser::getTaxaNames() const {
	return taxaNames;
}

void NewickParser::defineNamesRecursively() {
	root.memorizeNamesRecursively(taxaNames);
}

void NewickParser::defineInitialBranchLengthRecursively() {
	root.memorizeBLRecursively(initialBLs);
}

const std::map<size_t, double>& NewickParser::getInitialBLs() const {
	return initialBLs;
}

} /* namespace NewickReader */
} /* namespace Phylogeny */
