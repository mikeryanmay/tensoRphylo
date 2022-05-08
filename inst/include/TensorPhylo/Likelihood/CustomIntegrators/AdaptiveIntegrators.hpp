/*
 * AdaptiveIntegrators.hpp
 *
 *  Created on: Nov 24, 2019
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_CUSTOMINTEGRATORS_ADAPTIVEINTEGRATORS_HPP_
#define LIKELIHOOD_CUSTOMINTEGRATORS_ADAPTIVEINTEGRATORS_HPP_

#include <vector>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <iostream>

#include "Utils/Output/OutputManager.h"

//#define ANALYZE_ODE

namespace boost {
namespace numeric {
namespace odeint {

namespace detail {

template< class Stepper , class System , class State , class Time >
size_t integrate_adaptive_custom(
        Stepper stepper , System system , State &start_state ,
        Time &start_time , Time end_time , Time &dt ,
		std::vector<Time> &vecDeltaT, std::vector<Time> &vecTimes,  controlled_stepper_tag
)
{
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;

    //bool smallIntervalTried = false;
    failed_step_checker fail_checker;  // to throw a runtime_error if step size adjustment fails
    size_t count = 0;
    while( less_with_sign( start_time , end_time , dt ) ) {
        if( less_with_sign( end_time , static_cast<Time>(start_time + dt) , dt ) ) {
            dt = end_time - start_time;
        }
        if(dt < 0 && (start_time+dt) < 0) {
        	assert(fabs((double)(start_time+dt)) < std::numeric_limits<Time>::epsilon());
        	dt = end_time - start_time;
        }

        /*if((end_time - (start_time + dt) < 1.e-2) && !smallIntervalTried) {
        	std::cout << start_time << " -- " << end_time << " -- " << dt << std::endl;
        	dt = end_time - start_time;
        	smallIntervalTried = true;
        }*/

        controlled_step_result res;
        do {
        	//std::cout << "Trying : " << dt;
        	double attemptedDT = dt;
            res = st.try_step( system , start_state , start_time , dt );
            if(res == success) {
            	//std::cout << " -- Success : " << dt << std::endl;
            	start_state.roundNegativeProbabilityToZero();
            	vecDeltaT.push_back(attemptedDT);
            	// Start is incremented after a succes (start_time=start_time+dt)
            	vecTimes.push_back(start_time);

#ifdef ANALYZE_ODE
				std::ofstream oFile("tmpDmp.txt", std::ios::app);
				Likelihood::Monitor::ProbeState probe;
				probe = start_state.toProbe();
				probe.time = start_time + dt;
				oFile << probe.toDumpFormat();
#endif

            } else {
            	//std::cout << " -- Failure : " << dt << std::endl;
				vecDeltaT.push_back(attemptedDT);
            }
            fail_checker();  // check number of failed steps
        } while( res == fail );
        fail_checker.reset();  // if we reach here, the step was successful -> reset fail checker

        ++count;
    }
    return count;
}


template< class Stepper , class System , class State , class Time >
double integrate_adaptive_stepwise(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt ,
        bool doReset, dense_output_stepper_tag )
{
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;

    if(doReset) {
    	st.initialize( start_state , start_time , dt );
    }

    if( less_with_sign( st.current_time() , end_time , st.current_time_step() ) )
    {
    	bool isEndpoint = (st.current_time_step() > 0 && static_cast<Time>(st.current_time() + st.current_time_step()) >= end_time) ||
    			(st.current_time_step() < 0 && static_cast<Time>(st.current_time() + st.current_time_step()) <= end_time);
        if( !isEndpoint )
        {   //make sure we don't go beyond the end_time
            st.do_step( system );
        } else {
            // calculate time step to arrive exactly at end time
        	st.initialize( st.current_state() , st.current_time() , static_cast<Time>(end_time - st.current_time()) );
            st.do_step( system );
        }
    	start_state.roundNegativeProbabilityToZero();
    }
    // overwrite start_state with the final point
    boost::numeric::odeint::copy( st.current_state() , start_state );
    return st.current_time();
}


}

template< class Stepper , class System , class State , class Time >
size_t integrate_adaptive_custom(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt, std::vector<Time> &vecDeltaT, std::vector<Time> &vecTimes )
{
    typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;
    return detail::integrate_adaptive_custom(
            stepper , system , start_state ,
            start_time , end_time , dt ,
			vecDeltaT , vecTimes, stepper_category() );
}

template< class Stepper , class System , class State , class Time >
double integrate_adaptive_stepwise(
        Stepper stepper , System system , State &start_state ,
        Time start_time , Time end_time , Time dt, bool doReset )
{
    typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;
    return detail::integrate_adaptive_stepwise(
            stepper , system , start_state ,
            start_time , end_time , dt ,
			doReset, stepper_category() );
}

}
}
}


#endif /* LIKELIHOOD_CUSTOMINTEGRATORS_ADAPTIVEINTEGRATORS_HPP_ */
