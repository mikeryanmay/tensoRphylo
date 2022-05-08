/**
 * @file NewickGrammar.h
 *
 * @date Aug 22 2019
 * @author meyerx
 * @brief
 */
#ifndef NEWICKGRAMMAR_H_
#define NEWICKGRAMMAR_H_

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>

#include <string>

#include "TreeNode.h"

namespace Phylogeny {
namespace NewickReader {

namespace qi = ::boost::spirit::qi;

template<typename Iterator>
class NewickGrammar: public qi::grammar<Iterator, TreeNode()> {
public:
	NewickGrammar() : NewickGrammar::base_type(tree) {

		// For label use %= to assign the result of the parse to the string
		label %= qi::lexeme[+(~qi::char_(';') - ':' - ')' - ',' - '[')];

        // For branch length use %= to assign the result of the parse to the double
		branch_length %= ':' >> qi::double_;

		// Metadata for nodes
		metadata %= '&' >> qi::lexeme[+(~qi::char_(']') - '&')];

	    // When parsing the subtree just assign the elements that have been built in the subrules
		subtree = -descendant_list >> -label >> -( '[' >> +( metadata ) >> ']' ) >> -branch_length;

        // Descendant list is a vector of TreeNode, we just push back the created TreeNode into the vector
		descendant_list = '(' >> subtree >> *(',' >> subtree) >> ')' ;

		 // The tree receive the whole subtree using %=
		tree %= subtree >> ';';

		BOOST_SPIRIT_DEBUG_NODE(label);
		BOOST_SPIRIT_DEBUG_NODE(metadata);
		BOOST_SPIRIT_DEBUG_NODE(branch_length);
		BOOST_SPIRIT_DEBUG_NODE(subtree);
		BOOST_SPIRIT_DEBUG_NODE(descendant_list);
		BOOST_SPIRIT_DEBUG_NODE(tree);
	}

private:

	/* grammar rules */
	qi::rule<Iterator, TreeNode()> tree, subtree;
	qi::rule<Iterator, TreeNode_children()> descendant_list;
	qi::rule<Iterator, double()> branch_length;
	qi::rule<Iterator, std::string()> label;
	qi::rule<Iterator, std::vector<std::string>()> metadata;
};

} /* namespace NewickReader */
} /* namespace Phylogeny */

#endif /* NEWICKGRAMMAR_H_ */
