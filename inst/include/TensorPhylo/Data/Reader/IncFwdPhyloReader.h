/*
 * IncPhyloReader.h
 *
 *  Created on: Aug 28, 2019
 *      Author: xaviermeyer
 */

#ifndef DATA_READER_INCFWDPHYLOREADER_H_
#define DATA_READER_INCFWDPHYLOREADER_H_

#include <boost/smart_ptr/shared_ptr.hpp>

class NxsSimpleNode;

namespace Phylogeny {

namespace Data {
	class Container;
	typedef boost::shared_ptr<Container> ContainerSharedPtr;
}

namespace NexusReader {
	class NexusParser;
	typedef boost::shared_ptr<NexusParser> NexusParserSharedPtr;
}

namespace NewickReader {
	class NewickParser;
	typedef boost::shared_ptr<NewickParser> NewickParserSharedPtr;

	class TreeNode;

}

}

#endif /* DATA_READER_INCPHYLOREADER_H_ */
