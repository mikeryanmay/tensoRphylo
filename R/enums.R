#' Conditional probabilities
#'
#' @description
#' A list of conditional probability options for tensorphylo
#'
#' @details
#' These variables to specify the type of conditional probability calculations.
#' There are nine options:
#' | **Setting** | **Index** | **Interpretation** |
#' |:--------|:-:|:--------------|
#' | `conditionalProbability$TIME` _(default)_       | 0 | Don't condition on anything.
#' | `conditionalProbability$ROOT_SURVIVAL`          | 1 | Assume that the two lineages at the root each had to leave at least one extant sample.
#' | `conditionalProbability$ROOT_MRCA`              | 2 | Assume that the two lineages at the root each had to leave at least one extant sample, _AND_ that there was a speciation event at the root. **This is equivalent to conditioning on survival in most other programs.**
#' | `conditionalProbability$ROOT_SAMPLING`          | 3 | Assume that the two lineages at the root each had to leave at least one sample (extant or extinct).
#' | `conditionalProbability$ROOT_SAMPLING_AND_MRCA` | 8 | Assume that the two lineages at the root each had to leave at least one sample (extant or extinct), _AND_ that there was a speciation event at the root.
#' | `conditionalProbability$STEM_SURVIVAL`          | 4 | Assume that the one lineage at the stem had to leave at least one extant sample This is equivalent to conditioning on survival (with a stem) in most programs.
#' | `conditionalProbability$STEM_ONE_SAMPLE`        | 5 | Assume that the one lineage at the stem had to leave at least one sample (extant or extinct).
#' | `conditionalProbability$STEM_TWO_EXT_SAMPLES`   | 6 | Assume that the one lineage at the stem had to leave at least two extant samples.
#' | `conditionalProbability$STEM_TWO_SAMPLES`       | 7 | Assume that the one lineage at the stem had to leave at least two samples (extant or extinct).
#'
#' @examples
#' # create an empty TensorPhyloInstance object
#' tp <- new(TensorPhyloInstance, 4)
#'
#' # condition on the root speciation event and survival of the two descendants
#' tp$setConditionalProbabilityType(conditionalProbability$ROOT_MRCA)
#'
#' # this is equivalent to:
#' tp$setConditionalProbabilityType(1)
#' @name conditionalProbability
#' @export
NULL

conditionalProbability <- list(
  TIME=0,
	ROOT_SURVIVAL=1,
	ROOT_MRCA=2,
	ROOT_SAMPLING=3,
	STEM_SURVIVAL=4,
	STEM_ONE_SAMPLE=5,
	STEM_TWO_EXT_SAMPLES=6,
	STEM_TWO_SAMPLES=7,
	ROOT_SAMPLING_AND_MRCA=8
)

#' Approximation algorithms
#'
#' @description
#' A list of approximators for tensorphylo
#'
#' @details
#' These variables to specify the type of approximator to use.
#' There are five options:
#' | **Setting** | **Index** | **Interpretation** |
#' |:--------|:-:|:--------------|
#' | `approximatorVersion$AUTO_TUNING` _(default)_  | 0 | Use the approximator that is predicted to work the best for given tree size/state space/number of processors. Based on performance in simulation.
#' | `approximatorVersion$SEQUENTIAL_OPTIMIZED`     | 1 | Integrate all contemporaneous probabilities at each time point (i.e., concatenate linear operations across contemporaneous branches) (no parallelism).
#' | `approximatorVersion$SEQUENTIAL_BRANCHWISE`    | 2 | Integrate each branch independently (i.e., branch-specific probabilities are not integrated at the same time points) (no parallelism).
#' | `approximatorVersion$PARALLEL_OPTIMIZED`       | 3 | Integrate all contemporaneous probabilities at each time point (i.e., concatenate linear operations across contemporaneous branches) (linear operations are parallelized).
#' | `approximatorVersion$PARALLEL_BRANCHWISE`      | 4 | Integrate each branch independently (i.e., branch-specific probabilities are not integrated at the same time points) (linear operations are parallelized).
#'
#' @examples
#' # create an empty TensorPhyloInstance object
#' tp <- new(TensorPhyloInstance, 4)
#'
#' # specify a debug setting
#' tp$setLikelihoodApproximator(approximatorVersion$AUTO_TUNING)
#'
#' # this is equivalent to:
#' tp$setLikelihoodApproximator(0)
#' @name approximatorVersion
#' @export
NULL

approximatorVersion <- list(
  AUTO_TUNING=0,
  SEQUENTIAL_OPTIMIZED=1,
  SEQUENTIAL_BRANCHWISE=2,
  PARALLEL_OPTIMIZED=3,
  PARALLEL_BRANCHWISE=4
)

#' Integration algorithms
#'
#' @description
#' A list of numerical integration schemes for tensorphylo
#'
#' @details
#' These variables to specify the type of numerical integration algorithm to use.
#' There are four options:
#' | **Setting** | **Index** | **Interpretation** |
#' |:--------|:-:|:--------------|
#' | `integrationScheme$EULER`                      | 0 | The explicit Euler integrator.
#' | `integrationScheme$RUNGE_KUTTA4`               | 1 | The Runge-Kutta 4 integrator.
#' | `integrationScheme$RUNGE_KUTTA54` _(default)_  | 2 | The Runge-Kutta 4(5) integrator (AKA the Runge-Kutta-Fehlberg algorithm).
#' | `integrationScheme$RUNGE_KUTTA_DOPRI5`         | 3 | The Dormand–Prince integrator.
#'
#' @examples
#' # create an empty TensorPhyloInstance object
#' tp <- new(TensorPhyloInstance, 4)
#'
#' # specify a debug setting
#' tp$setIntegrationScheme(integrationScheme$RUNGE_KUTTA_DOPRI5)
#'
#' # this is equivalent to:
#' tp$setIntegrationScheme(3)
#' @name integrationScheme
#' @export
NULL

integrationScheme <- list(
	EULER=0,
	RUNGE_KUTTA4=1,
	RUNGE_KUTTA54=2,
	RUNGE_KUTTA_DOPRI5=3
)

#' TensorPhyloInstance debug modes
#'
#' @description
#' A list of debug modes for tensorphylo
#'
#' @details
#' These variables to specify debug modes
#' There are three options:
#' | **Setting** | **Index** | **Interpretation** |
#' |:--------|:-:|:--------------|
#' | `debugMode$DBG_NONE` _(default)_  | 0 | No debug info
#' | `debugMode$DBG_PRINT`             | 1 | Debug info printed to screen
#' | `debugMode$DBG_FILE`              | 2 | Debug info printed to file
#'
#' @examples
#' # create an empty TensorPhyloInstance object
#' tp <- new(TensorPhyloInstance, 4)
#'
#' # specify a debug setting
#' tp$setDebugMode(debugMode$DBG_PRINT)
#'
#' # this is equivalent to:
#' tp$setDebugMode(1)
#' @name debugMode
#' @export
NULL

debugMode <- list(
  DBG_NONE=0,
  DBG_PRINT=1,
  DBG_FILE=2
)
