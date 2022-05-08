/*
 * AsyncParameterContainer.cpp
 *
 *  Created on: March 10, 2020
 *      Author: xaviermeyer
 */

#include "AsyncParameterContainer.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>

namespace Parameters {

AsyncParameterContainer::AsyncParameterContainer() : nState(0) {
}

AsyncParameterContainer::~AsyncParameterContainer() {
}

} /* namespace Parameters */
