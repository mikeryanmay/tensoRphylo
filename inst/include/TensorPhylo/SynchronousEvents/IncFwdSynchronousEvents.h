/*
 * IncFwdSynchronousEvent.h
 *
 *  Created on: Sep 3, 2019
 *      Author: xaviermeyer
 */

#ifndef SYNCHRONOUSEVENT_INCFWDSYNCHRONOUSEVENT_H_
#define SYNCHRONOUSEVENT_INCFWDSYNCHRONOUSEVENT_H_

#include <boost/smart_ptr/shared_ptr.hpp>

namespace SynchronousEvents {

	class Definition;
	typedef boost::shared_ptr<Definition> DefinitionSharedPtr;

	class Container;
	typedef boost::shared_ptr<Container> ContainerSharedPtr;

}



#endif /* SYNCHRONOUSEVENT_INCFWDSYNCHRONOUSEVENT_H_ */
