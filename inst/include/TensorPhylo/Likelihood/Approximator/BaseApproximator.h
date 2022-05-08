/*
 * BaseApproximator.h
 *
 *  Created on: Aug 25, 2019
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_APPROXIMATOR_BASEAPPROXIMATE_H_
#define LIKELIHOOD_APPROXIMATOR_BASEAPPROXIMATE_H_


#include <vector>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "Tensor/IncFwdTensor.h"
#include "SynchronousEvents/IncFwdSynchronousEvents.h"
#include "Data/Reader/IncFwdPhyloReader.h"
#include "Likelihood/Monitor/ProbeState.h"
#include "Likelihood/Scheduler/IncFwdScheduler.h"
#include "Likelihood/ConditionTypes/ConditionType.h"
#include "Likelihood/CustomIntegrators/IntegratorFactory.hpp"

namespace Likelihood {
namespace Approximator {

class BaseApproximator {
public:
	BaseApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
					 Conditions::conditionalProbability_t aConditionType,
					 Phylogeny::Data::ContainerSharedPtr aPtrData,
					 Scheduler::SchedulerSharedPtr aPtrScheduler,
					 SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
					 Tensor::ContainerSharedPtr aPtrTensorCont);
	virtual ~BaseApproximator();

	virtual void enableTreeLikelihoodCorrection();
	virtual void disableTreeLikelihoodCorrection();

	// also implement these in AutoTuningApproximator
	virtual void setQuasistationaryFrequencyMode(bool setActive);
	virtual Eigen::VectorXd getRootFrequency(double t);
	virtual void setPriorStateProbability(const Eigen::VectorXd &aPriorStateProbability);

	virtual void setConditionalProbabilityCompatibilityMode(bool setActive);

	double approximateLikelihood();
	virtual double approximateLogLikelihood() = 0;

	virtual void orderProbes();
	virtual const std::vector<Likelihood::Monitor::ProbeState>& getObservedProbesState() const;
	virtual const std::vector<double>& getIntegrationTimes() const;

	virtual void setDefaultDeltaT(double aDeltaT) = 0;
	virtual size_t getTotalNumberOfIntegrationSteps() const = 0;

public:
	static const double DEFAULT_DELTA_T, DEFAULT_ABS_TOLERANCE, DEFAULT_REL_TOLERANCE;

protected:

	bool applyTreeCorrection, useQuasistationaryFrequency, condCompatibilityMode;
	double deltaT, logLikelihood;
	boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::mean > > deltaTAcc;

	Likelihood::Integrator::integrationScheme_t intScheme;
	Likelihood::Conditions::conditionalProbability_t conditionType;
	Eigen::VectorXd priorStateProbability;

	Phylogeny::Data::ContainerSharedPtr  ptrData;
	Scheduler::SchedulerSharedPtr ptrScheduler;
	SynchronousEvents::ContainerSharedPtr ptrSyncEventsCont;
	Tensor::ContainerSharedPtr ptrTensorCont;

	std::vector<double> integrationTimes;
	std::vector<Likelihood::Monitor::ProbeState> vecProbesState;

	bool areEventsPossible();

};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_APPROXIMATOR_BASEAPPROXIMATE_H_ */
