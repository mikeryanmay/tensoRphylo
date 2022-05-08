/*
 * IntegratorFactory.h
 *
 *  Created on: Dec 6, 2019
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_CUSTOMINTEGRATORS_INTEGRATORFACTORY_H_
#define LIKELIHOOD_CUSTOMINTEGRATORS_INTEGRATORFACTORY_H_
#include <cstddef>
#include <vector>

#include "CustomDenseRK.hpp"

#include <boost/numeric/odeint.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/rolling_count.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>

namespace Likelihood {
namespace Integrator {

typedef enum {
	EULER=0,
	RUNGE_KUTTA4=1,
	RUNGE_KUTTA54=2,
	RUNGE_KUTTA_DOPRI5=3,
	DENSE_RUNGE_KUTTA_DOPRI5=4
} integrationScheme_t;

/************************************************/
/*************** BASE INTEGRATOR ****************/
/************************************************/
template <class StateType, class IntegratorKernel, class OperationType>
class Base {
public:
	Base(const double aAbsError, const double aRelError, const double aDeltaT);
	virtual ~Base();

	void setDeltaT(double aDeltaT);
	size_t getNSteps() const;

	virtual int integrate(double startTime, double endTime,
						  StateType &state, IntegratorKernel &intKernel) = 0;
	virtual void reset() = 0;

	const std::vector<double>& getVecTimes() const;

protected:
	static const size_t DELTA_T_ACCUMULATOR_WINDOW_SIZE;
	double DELTA_T;
	const double ABS_ERROR, REL_ERROR;
	size_t nSteps;
	std::vector<double> vecTimes;

};

/************************************************/
/*****************    EULER    ******************/
/************************************************/

template <class StateType, class IntegratorKernel, class OperationType>
class Euler: public Base<StateType, IntegratorKernel, OperationType> {
public:
	Euler(const double aAbsError, const double aRelError, const double aDeltaT);
	~Euler();

	int integrate(double startTime, double endTime,
				  StateType &state, IntegratorKernel &intKernel);
	void reset();

private:
	boost::numeric::odeint::euler< StateType ,
								   double ,
								   StateType,
								   double ,
								   boost::numeric::odeint::vector_space_algebra,
								   OperationType,
								   boost::numeric::odeint::always_resizer > constStepper;
};

/************************************************/
/*************    Runge Kutta 4    **************/
/************************************************/

template <class StateType, class IntegratorKernel, class OperationType>
class RungeKutta4: public Base<StateType, IntegratorKernel, OperationType> {
public:
	RungeKutta4(const double aAbsError, const double aRelError, const double aDeltaT);
	~RungeKutta4();

	int integrate(double startTime, double endTime,
				  StateType &state, IntegratorKernel &intKernel);
	void reset();

private:
	boost::numeric::odeint::runge_kutta4< StateType ,
										  double ,
										  StateType,
										  double ,
										  boost::numeric::odeint::vector_space_algebra,
										  OperationType,
										  boost::numeric::odeint::always_resizer > constStepper;
};

/************************************************/
/************    Runge Kutta 54    **************/
/************************************************/

template <class StateType, class IntegratorKernel, class OperationType>
class RungeKutta54: public Base<StateType, IntegratorKernel, OperationType> {
public:
	RungeKutta54(const double aAbsError, const double aRelError, const double aDeltaT);
	~RungeKutta54();

	int integrate(double startTime, double endTime,
				  StateType &state, IntegratorKernel &intKernel);
	void reset();

private:
	typedef boost::numeric::odeint::runge_kutta_cash_karp54< StateType ,
															 double ,
															 StateType,
															 double ,
															 boost::numeric::odeint::vector_space_algebra,
															 OperationType,
															 boost::numeric::odeint::always_resizer > rk54_stepper_t;
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::rolling_count, boost::accumulators::tag::rolling_mean > > deltaTAcc;
	typedef boost::numeric::odeint::controlled_runge_kutta<rk54_stepper_t> controlledRK54_t;
	controlledRK54_t adaptiveStepper;
};

/************************************************/
/*********    Runge Kutta DOPRI 5    ************/
/************************************************/

template <class StateType, class IntegratorKernel, class OperationType>
class RungeKuttaDOPRI5: public Base<StateType, IntegratorKernel, OperationType> {
public:
	RungeKuttaDOPRI5(const double aAbsError, const double aRelError, const double aDeltaT);
	~RungeKuttaDOPRI5();

	int integrate(double startTime, double endTime,
				  StateType &state, IntegratorKernel &intKernel);
	void reset();
private:
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::rolling_count, boost::accumulators::tag::rolling_mean > > deltaTAcc;
	typedef boost::numeric::odeint::runge_kutta_dopri5< StateType ,
													 	double ,
														StateType,
														double ,
														boost::numeric::odeint::vector_space_algebra,
														OperationType,
														boost::numeric::odeint::always_resizer > rkd5_stepper_t;
	typedef boost::numeric::odeint::controlled_runge_kutta<rkd5_stepper_t> controlledRDK5_t;
	controlledRDK5_t adaptiveStepper;
};

/************************************************/
/*********    Runge Kutta DOPRI 5    ************/
/************************************************/

template <class StateType, class IntegratorKernel, class OperationType>
class DenseRungeKuttaDOPRI5: public Base<StateType, IntegratorKernel, OperationType> {
public:
	typedef boost::numeric::odeint::runge_kutta_dopri5< StateType ,
													 	double ,
														StateType,
														double ,
														boost::numeric::odeint::vector_space_algebra,
														OperationType,
														boost::numeric::odeint::always_resizer > rkd5_stepper_t;
	typedef boost::numeric::odeint::controlled_runge_kutta<rkd5_stepper_t> controlledRDK5_t;

	typedef boost::numeric::odeint::custom_dense_runge_kutta< controlledRDK5_t> customDenseRDK5_t;

public:
	DenseRungeKuttaDOPRI5(const double aAbsError, const double aRelError, const double aDeltaT);
	~DenseRungeKuttaDOPRI5();

	int integrate(double startTime, double endTime,
				  StateType &state, IntegratorKernel &intKernel);
	void reset();

	const std::vector<rkd5_stepper_t>& getSteppers() const;
	const std::vector< std::pair<StateType, StateType> >& getStepStates() const;
	const std::vector< std::pair<StateType, StateType> >& getStepDerivs() const;
	const std::vector< std::pair<double, double> >& getStepTimes() const;

	void transferStepsMemory(std::vector< rkd5_stepper_t >& aSteppers,
							 std::vector< std::pair< double, double > >& aTimes,
							 std::vector< std::pair< StateType, StateType >  >& aStepStates,
							 std::vector< std::pair< StateType, StateType >  >& aStepDerivs );

private:

	customDenseRDK5_t adaptiveStepper;

};


/************************************************/
/***************    FACTORY    ******************/
/************************************************/

namespace Factory {

template <class StateType, class IntegratorKernel, class OperationType>
Base<StateType, IntegratorKernel, OperationType>* createIntegrator(const double aAbsError, const double aRelError, const double aDeltaT, integrationScheme_t aIntScheme) {
	if(aIntScheme == EULER) {
		return new Euler<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT);
	} else if(aIntScheme == RUNGE_KUTTA4) {
		return new RungeKutta4<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT);
	} else if(aIntScheme == RUNGE_KUTTA54) {
		return new RungeKutta54<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT);
	} else if(aIntScheme == RUNGE_KUTTA_DOPRI5) {
		return new RungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT);
	} else if(aIntScheme == DENSE_RUNGE_KUTTA_DOPRI5) {
		return new DenseRungeKuttaDOPRI5<StateType, IntegratorKernel, OperationType>(aAbsError, aRelError, aDeltaT);
	} else {
		assert(false && "Integrator unknown from the Integrator factory.");
		return NULL;
	}
}

}

} /* namespace Integrator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_CUSTOMINTEGRATORS_INTEGRATORFACTORY_H_ */
