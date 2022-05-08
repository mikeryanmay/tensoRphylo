/*
 * AutoTuningApproximator.h
 *
 *  Created on: Apr 21, 2020
 *      Author: meyerx
 */

#ifndef LIKELIHOOD_AUTOTUNINGAPPROXIMATOR_H_
#define LIKELIHOOD_AUTOTUNINGAPPROXIMATOR_H_

#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>
#include "BaseApproximator.h"

#include "IncLikelihoodApproximator.h"
#include "Utils/Profiling/CustomProfiling.h"


#include <map>
#include <utility>
#include <vector>


namespace Likelihood {
namespace Approximator {

class AutoTuningApproximator : public BaseApproximator {
public:
	AutoTuningApproximator(Likelihood::Integrator::integrationScheme_t aIntScheme,
   	   	   	Conditions::conditionalProbability_t aConditionType,
			Phylogeny::Data::ContainerSharedPtr aPtrData,
			Scheduler::SchedulerSharedPtr aPtrScheduler,
			SynchronousEvents::ContainerSharedPtr aPtrSyncEventsCont,
			Tensor::ContainerSharedPtr aPtrTensorCont);
	~AutoTuningApproximator();

	void enableTreeLikelihoodCorrection();
	void disableTreeLikelihoodCorrection();

	void setQuasistationaryFrequencyMode(bool setActive);
	Eigen::VectorXd getRootFrequency(double t);
	void setPriorStateProbability(const Eigen::VectorXd &aPriorStateProbability);

	void setConditionalProbabilityCompatibilityMode(bool setActive);

	void orderProbes();
	const std::vector<Likelihood::Monitor::ProbeState>& getObservedProbesState() const;
	const std::vector<double>& getIntegrationTimes() const;

	double approximateLogLikelihood();

	void setDefaultDeltaT(double aDeltaT);
	size_t getTotalNumberOfIntegrationSteps() const;

private:

	static const size_t ROLLING_MEAN_WINDOW_SIZE, AUTO_TUNING_LENGTH;

	approximatorVersion_t currentApproximatorType;
	ApproximatorSharedPtr ptrApprox;

	double avgTimeAtSelection;
	size_t likCounter;

	Utils::Profiling::CustomProfiling cp;
	typedef boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::rolling_count, boost::accumulators::tag::rolling_mean > > accumulator_t;
	std::vector<accumulator_t> nThreadPerf;

	static const size_t thresholdNTips, thresholdNStates;
	static const std::vector<double> coeffsSmall1A_1T, coeffsSmall2A_1T, coeffsSmall3A_2T, coeffsSmall3A_4T, coeffsSmall4A_2T, coeffsSmall4A_4T;
	static const std::vector<double> coeffsLarge1A_1T, coeffsLarge2A_1T, coeffsLarge3A_2T, coeffsLarge3A_4T, coeffsLarge4A_2T, coeffsLarge4A_4T;
	static const std::map< std::pair<size_t, size_t>, std::vector<double> > coeffsSmall, coeffsLarge;

	static const size_t thresholdNTips_clado, thresholdNStates_clado;
	static const std::vector<double> coeffsSmall2A_1T_clado, coeffsSmall1A_1T_clado, coeffsSmall4A_2T_clado, coeffsSmall3A_2T_clado, coeffsSmall4A_4T_clado, coeffsSmall3A_4T_clado;
	static const std::vector<double> coeffsLarge2A_1T_clado, coeffsLarge1A_1T_clado, coeffsLarge4A_2T_clado, coeffsLarge3A_2T_clado, coeffsLarge4A_4T_clado, coeffsLarge3A_4T_clado;
	static const std::map< std::pair<size_t, size_t>, std::vector<double> > coeffsSmall_clado, coeffsLarge_clado;


	static const size_t thresholdNTips_quasse, thresholdNStates_quasse;
	static const std::vector<double> coeffsSmall1A_1T_quasse, coeffsSmall2A_1T_quasse, coeffsSmall3A_2T_quasse, coeffsSmall4A_2T_quasse;
	static const std::vector<double> coeffsLarge1A_1T_quasse, coeffsLarge2A_1T_quasse, coeffsLarge3A_2T_quasse, coeffsLarge4A_2T_quasse;
	static const std::map< std::pair<size_t, size_t>, std::vector<double> > coeffsSmall_quasse, coeffsLarge_quasse;


	bool isLikelihoodDroppingWithNT(size_t currentNThreads, size_t maxNThreads);
	size_t defineNextSettingToTry(size_t currentNThreads, size_t maxNThreads);


	double predictLikelihoodEvaluationTimeMusse(size_t nTips, size_t kStates, const std::vector<double> &coeffs) const;
	double predictLikelihoodEvaluationTimeClasse(size_t nTips, size_t kStates, const std::vector<double> &coeffs) const;
	double predictLikelihoodEvaluationTimeQuasse(size_t nTips, size_t kStates, const std::vector<double> &coeffs) const;
	approximatorVersion_t predictBestIntegrator();

	ApproximatorSharedPtr createApproximator(approximatorVersion_t anApproximatorType);


};

} /* namespace Approximator */
} /* namespace Likelihood */

#endif /* LIKELIHOOD_AUTOTUNINGAPPROXIMATOR_H_ */
