/*
 * InclFwdTreeStructure.h
 *
 *  Created on: Aug 25, 2019
 *      Author: meyerx
 */

#ifndef DATA_STRUCTURE_INCFWDTREESTRUCTURE_H_
#define DATA_STRUCTURE_INCFWDTREESTRUCTURE_H_

#include <boost/smart_ptr/shared_ptr.hpp>

namespace Phylogeny {



namespace Structure {
	class Edge;
	class Node;
	class Tree;
	typedef boost::shared_ptr<Tree> TreeSharedPtr;

}
}

namespace PS = ::Phylogeny::Structure;


#endif /* DATA_STRUCTURE_INCFWDTREESTRUCTURE_H_ */
