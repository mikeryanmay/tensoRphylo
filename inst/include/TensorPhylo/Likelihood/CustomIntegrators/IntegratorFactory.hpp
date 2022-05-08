/*
 * IntegratorFactory.hpp
 *
 *  Created on: Dec 6, 2019
 *      Author: meyerx
 */

#include "IntegratorFactory.h"

#include "AdaptiveIntegrators.hpp"

#include <boost/numeric/odeint.hpp>

#ifndef LIKELIHOOD_CUSTOMINTEGRATORS_INTEGRATORFACTORY_DEF_H_
#define LIKELIHOOD_CUSTOMINTEGRATORS_INTEGRATORFACTORY_DEF_H_

namespace Likelihood {
namespace Integrator {

/************************************************/
/*************** BASE INTEGRATOR ****************/
/************************************************/

template <class StateType, class IntegratorKernel, class OperationType>
const size_t Base<StateType, IntegratorKernel, OperationType>::DELTA_T_ACCUMULATOR_WINDOW_SIZE = 100;

template <class StateType, class IntegratorKernel, class OperationType>
Base<StateType, IntegratorKernel, OperationType>::Base(const double aAbsError, const double aRelError, const double aDeltaT) :
												  DELTA_T(aDeltaT), ABS_ERROR(aAbsError), REL_ERROR(aRelError) {
	nSteps = 0;
}

template <class StateType, class IntegratorKernel, class OperationType>
Base<StateType, IntegratorKernel, OperationType>::~Base() {
}

template <class StateType, class IntegratorKernel, class OperationType>
void Base<StateType, IntegratorKernel, OperationType>::setDeltaT(double aDeltaT) {
	DELTA_T = aDeltaT;
}

template <class StateType, class IntegratorKernel, class OperationType>
size_t Base<StateType, IntegratorKernel, OperationType>::getNSteps() const {
	return nSteps;
}

template <class StateType, class IntegratorKernel, class OperationType>
const std::vector<double>& Base<StateType, IntegratorKernel, OperationType>::getVecTimes() const {
	return vecTimes;
}


/************************************************/
/*****************    EULER    ******************/
/************************************************/
template <class StateType, class IntegratorKernel, class OperationType>
Euler<StateType, IntegratorKernel, OperationType>::Euler(const double aAbsError, const double aRelError, const double aDeltaT) :
												   Base<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT) {
}

template <class StateType, class IntegratorKernel, class OperationType>
Euler<StateType, IntegratorKernel, OperationType>::~Euler() {
}

template <class StateType, class IntegratorKernel, class OperationType>
int Euler<StateType, IntegratorKernel, OperationType>::integrate(double startTime, double endTime,
																 StateType &state, IntegratorKernel &intKernel) {
	int steps = integrate_const( constStepper, boost::ref(intKernel), state, startTime, endTime, Base<StateType, IntegratorKernel, OperationType>::DELTA_T);
	state.roundNegativeProbabilityToZero();
	Base<StateType, IntegratorKernel, OperationType>::nSteps += steps;
	return steps;
}

template <class StateType, class IntegratorKernel, class OperationType>
void Euler<StateType, IntegratorKernel, OperationType>::reset() {
	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
}

/************************************************/
/*************    Runge Kutta 4    **************/
/************************************************/
template <class StateType, class IntegratorKernel, class OperationType>
RungeKutta4<StateType, IntegratorKernel, OperationType>::RungeKutta4(const double aAbsError, const double aRelError, const double aDeltaT) :
														 Base<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT) {
}

template <class StateType, class IntegratorKernel, class OperationType>
RungeKutta4<StateType, IntegratorKernel, OperationType>::~RungeKutta4() {
}

template <class StateType, class IntegratorKernel, class OperationType>
int RungeKutta4<StateType, IntegratorKernel, OperationType>::integrate(double startTime, double endTime,
																	   StateType &state, IntegratorKernel &intKernel) {
	int steps = integrate_const( constStepper, boost::ref(intKernel), state, startTime, endTime, Base<StateType, IntegratorKernel, OperationType>::DELTA_T);
	state.roundNegativeProbabilityToZero();
	Base<StateType, IntegratorKernel, OperationType>::nSteps += steps;
	return steps;
}

template <class StateType, class IntegratorKernel, class OperationType>
void RungeKutta4<StateType, IntegratorKernel, OperationType>::reset() {
	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
}

/************************************************/
/************    Runge Kutta 45    **************/
/************************************************/
template <class StateType, class IntegratorKernel, class OperationType>
RungeKutta54<StateType, IntegratorKernel, OperationType>::RungeKutta54(const double aAbsError, const double aRelError, const double aDeltaT) :
														  Base<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT),
														  deltaTAcc(boost::accumulators::tag::rolling_window::window_size = Base<StateType, IntegratorKernel, OperationType>::DELTA_T_ACCUMULATOR_WINDOW_SIZE),
														  adaptiveStepper(make_controlled( Base<StateType, IntegratorKernel, OperationType>::ABS_ERROR,
																  	  	  	  	  	  	   Base<StateType, IntegratorKernel, OperationType>::REL_ERROR, rk54_stepper_t())) {
}

template <class StateType, class IntegratorKernel, class OperationType>
RungeKutta54<StateType, IntegratorKernel, OperationType>::~RungeKutta54() {
}

template <class StateType, class IntegratorKernel, class OperationType>
int RungeKutta54<StateType, IntegratorKernel, OperationType>::integrate(double startTime, double endTime,
																	    StateType &state, IntegratorKernel &intKernel) {

	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
	double deltaT = Base<StateType, IntegratorKernel, OperationType>::DELTA_T;
	if(boost::accumulators::rolling_count(deltaTAcc) > Base<StateType, IntegratorKernel, OperationType>::DELTA_T_ACCUMULATOR_WINDOW_SIZE/2) {
		deltaT =  boost::accumulators::rolling_mean(deltaTAcc);
	}

	std::vector<double> vecDeltaT;
	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
	int steps = integrate_adaptive_custom( boost::ref(adaptiveStepper), boost::ref(intKernel), state, startTime, endTime, deltaT, vecDeltaT, Base<StateType, IntegratorKernel, OperationType>::vecTimes);
	Base<StateType, IntegratorKernel, OperationType>::nSteps += steps;
	RungeKutta54<StateType, IntegratorKernel, OperationType>::reset();

	for(size_t iD=1; iD<vecDeltaT.size()-1; ++iD) {
		deltaTAcc(vecDeltaT[iD]);
	}

	return steps;
}

template <class StateType, class IntegratorKernel, class OperationType>
void RungeKutta54<StateType, IntegratorKernel, OperationType>::reset() {
	// hard reset
	/*adaptiveStepper = make_controlled(Base<StateType, IntegratorKernel, OperationType>::ABS_ERROR,
	  	  	  	   Base<StateType, IntegratorKernel, OperationType>::REL_ERROR, rk54_stepper_t());*/
	//adaptiveStepper.reset();
}

/************************************************/
/*********    Runge Kutta DOPRI 5    ************/
/************************************************/
template <class StateType, class IntegratorKernel, class OperationType>
RungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::RungeKuttaDOPRI5(const double aAbsError, const double aRelError, const double aDeltaT) :
															  Base<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT),
															  deltaTAcc(boost::accumulators::tag::rolling_window::window_size = Base<StateType, IntegratorKernel, OperationType>::DELTA_T_ACCUMULATOR_WINDOW_SIZE),
															  adaptiveStepper(make_controlled(Base<StateType, IntegratorKernel, OperationType>::ABS_ERROR,
																	                          Base<StateType, IntegratorKernel, OperationType>::REL_ERROR, rkd5_stepper_t())) {
}

template <class StateType, class IntegratorKernel, class OperationType>
RungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::~RungeKuttaDOPRI5() {
}

template <class StateType, class IntegratorKernel, class OperationType>
int RungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::integrate(double startTime, double endTime,
																		    StateType &state, IntegratorKernel &intKernel) {

	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
	double deltaT = Base<StateType, IntegratorKernel, OperationType>::DELTA_T;
	if(boost::accumulators::rolling_count(deltaTAcc) > Base<StateType, IntegratorKernel, OperationType>::DELTA_T_ACCUMULATOR_WINDOW_SIZE/2) {
		deltaT = boost::accumulators::rolling_mean(deltaTAcc);
	}

	std::vector<double> vecDeltaT;
	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
	int steps = integrate_adaptive_custom( boost::ref(adaptiveStepper), boost::ref(intKernel), state, startTime, endTime , deltaT, vecDeltaT, Base<StateType, IntegratorKernel, OperationType>::vecTimes);
	reset();

	Base<StateType, IntegratorKernel, OperationType>::nSteps += steps;

	for(size_t iD=1; iD<vecDeltaT.size()-1; ++iD) {
		deltaTAcc(vecDeltaT[iD]);
	}
	return steps;
}

template <class StateType, class IntegratorKernel, class OperationType>
void RungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::reset() {
	// hard reset
	adaptiveStepper.reset();
	//adaptiveStepper = make_controlled(Base<StateType, IntegratorKernel, OperationType>::ABS_ERROR,
	//				  Base<StateType, IntegratorKernel, OperationType>::REL_ERROR, rkd5_stepper_t());
}

/************************************************/
/****** Dense Runge Kutta DOPRI 5    **********/
/************************************************/
template <class StateType, class IntegratorKernel, class OperationType>
DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::DenseRungeKuttaDOPRI5(const double aAbsError, const double aRelError, const double aDeltaT) :
															  Base<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT),
															  adaptiveStepper(make_controlled(Base<StateType, IntegratorKernel, OperationType>::ABS_ERROR,
																	                       Base<StateType, IntegratorKernel, OperationType>::REL_ERROR, rkd5_stepper_t())) {

}

template <class StateType, class IntegratorKernel, class OperationType>
DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::~DenseRungeKuttaDOPRI5() {
}

template <class StateType, class IntegratorKernel, class OperationType>
int DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::integrate(double startTime, double endTime,
																		    StateType &state, IntegratorKernel &intKernel) {

	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
	double deltaT = Base<StateType, IntegratorKernel, OperationType>::DELTA_T;

	std::vector<double> vecDeltaT;
	Base<StateType, IntegratorKernel, OperationType>::vecTimes.clear();
	boost::numeric::odeint::integrate_adaptive( boost::ref(adaptiveStepper), boost::ref(intKernel), state, startTime, endTime, deltaT);
	//reset();

	Base<StateType, IntegratorKernel, OperationType>::nSteps = adaptiveStepper.getSteppers().size();

	for(size_t iT=0; iT < adaptiveStepper.getStepTimes().size(); ++iT ) {
		Base<StateType, IntegratorKernel, OperationType>::vecTimes.push_back(adaptiveStepper.getStepTimes()[iT].second);
	}

	return Base<StateType, IntegratorKernel, OperationType>::vecTimes.size();
}

template <class StateType, class IntegratorKernel, class OperationType>
void DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::reset() {
	//adaptiveStepper.reset();
}

template <class StateType, class IntegratorKernel, class OperationType>
const std::vector< typename DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::rkd5_stepper_t >& DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::getSteppers() const {
	return adaptiveStepper.getSteppers();
}

template <class StateType, class IntegratorKernel, class OperationType>
const std::vector< std::pair<StateType, StateType> >& DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::getStepStates() const {
	return adaptiveStepper.getStepStates();
}

template <class StateType, class IntegratorKernel, class OperationType>
const std::vector< std::pair<StateType, StateType> >& DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::getStepDerivs() const {
	return adaptiveStepper.getStepDerivs();
}

template <class StateType, class IntegratorKernel, class OperationType>
const std::vector< std::pair<double, double> >& DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::getStepTimes() const {
	return adaptiveStepper.getStepTimes();
}

template <class StateType, class IntegratorKernel, class OperationType>
void DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>::transferStepsMemory(
						 std::vector< rkd5_stepper_t >& aSteppers,
						 std::vector< std::pair< double, double > >& aTimes,
						 std::vector< std::pair< StateType, StateType >  >& aStepStates,
						 std::vector< std::pair< StateType, StateType >  >& aStepDerivs ) {
	adaptiveStepper.transferStepsMemory(boost::ref(aSteppers), aTimes, aStepStates, aStepDerivs);
}

/************************************************/
/***************    FACTORY    ******************/
/************************************************/

} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_CUSTOMINTEGRATORS_INTEGRATORFACTORY_DEF_H_ */
