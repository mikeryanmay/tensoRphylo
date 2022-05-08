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
/*
 * *SimplexCuda* is a LP solver using the simplex method on GPU.
 * Copyright (C) <2010>  <Xavier Meyer>
 *
 *  Created on: Sep 11, 2010
 *      Author: meyerx
 */

#include "CustomProfiling.h"

#include <stdio.h>

namespace Utils {
namespace Profiling {

CustomProfiling::CustomProfiling()  {

}

CustomProfiling::~CustomProfiling(){
}


void CustomProfiling::startTime(){
	 startTime(0);
}

void CustomProfiling::endTime(){
	endTime(0);
}

double CustomProfiling::duration() const {
	return duration(0);
}


void CustomProfiling::startTime(int iClock){
	 gettimeofday(&start[iClock], NULL);
}

void CustomProfiling::endTime(int iClock){
	 gettimeofday(&end[iClock], NULL);
}

double CustomProfiling::duration(int iClock) const {
	long s = end[iClock].tv_sec - start[iClock].tv_sec;
	long us = end[iClock].tv_usec - start[iClock].tv_usec;
	return ((double)s + us*1e-6);
}


} // Profiling
} // Utils


