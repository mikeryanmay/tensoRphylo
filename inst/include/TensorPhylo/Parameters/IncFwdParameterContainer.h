/*
 * IncFwdParameterContainer.h
 *
 *  Created on: Mar 10, 2020
 *      Author: xaviermeyer
 */

#ifndef PARAMETERS_INCFWDPARAMETERCONTAINER_H_
#define PARAMETERS_INCFWDPARAMETERCONTAINER_H_

#include <boost/smart_ptr/shared_ptr.hpp>

namespace Parameters {
	class Container;
	typedef boost::shared_ptr<Container> ContainerSharedPtr;

	class AsyncParameterContainer;
	typedef boost::shared_ptr<AsyncParameterContainer> AsyncContainerSharedPtr;

	class SyncParameterContainer;
	typedef boost::shared_ptr<SyncParameterContainer> SyncContainerSharedPtr;

}




#endif /* PARAMETERS_INCFWDPARAMETERCONTAINER_H_ */
