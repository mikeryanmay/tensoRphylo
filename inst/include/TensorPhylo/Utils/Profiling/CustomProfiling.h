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

#ifndef CUSTOMPROFLING_H_
#define CUSTOMPROFLING_H_

#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <string>

#include <iosfwd>

#if PROFILING
#define _START_EVENT_CP(profiler, iClock) 			profiler.startTime(iClock);
#define _END_EVENT_CP(profiler, iClock)				profiler.endTime(iClock);
#define _DURATION_EVENT_CP(profiler, iClock)		profiler.duration(iClock);
#else
#define _START_EVENT_CP(profiler, iClock) 		void();
#define _END_EVENT_CP(profiler, iClock)			void();
#define _DURATION_EVENT_CP(profiler, iClock)	0.;
#endif

namespace Utils {
namespace Profiling {

class CustomProfiling {
public:

	CustomProfiling();
	~CustomProfiling();

	void startTime();
	void endTime();
	double duration() const;

	void startTime(int iClock);
	void endTime(int iClock);
	double duration(int iClock) const;

private:

	timeval start[10];
	timeval end[10];

};

} // Profiling
} // Utils

#endif /* CUSTOMPROFLING_H_ */
