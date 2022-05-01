#' Conditional probabilities
#'
#' @description
#' A list of conditional probability options for tensorphylo
#'
#' @details
#' Use these variables to specify the type of conditional probability calculations.
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
#' # create an empty tensorphylo object
#' tp <- new(TensorPhylo, 4)
#'
#' # condition on survival of the root lineages
#' tp$setConditionalProbabilityType( conditionalProbability$ROOT_SURVIVAL )
#'
#' # this is equivalent to:
#' tp$setConditionalProbabilityType( 1 )
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
#' Use these variables to specify the type of approximator to use.
#' There are five options:
#' | **Setting** | **Index** | **Interpretation** |
#' |:--------|:-:|:--------------|
#' | `approximatorVersion$AUTO_TUNING` _(default)_  | 0 | Use the approximator that works the best after the first few calculations.
#' | `approximatorVersion$SEQUENTIAL_OPTIMIZED`     | 1 | Integrate each branch one at a time, but share extinction probability calculations.
#' | `approximatorVersion$SEQUENTIAL_BRANCHWISE`    | 2 | Integrate each branch one at a time, with separate extinction probability calculations.
#' | `approximatorVersion$PARALLEL_OPTIMIZED`       | 3 | Parallelize integration over branches, but share extinction probability calclations.
#' | `approximatorVersion$PARALLEL_BRANCHWISE`      | 4 | Parallelize integration over branches, with separate extinction probability calculations.
#'
#' @examples
#' # create an empty tensorphylo object
#' tp <- new(TensorPhylo, 4)
#'
#' # specify a debug setting
#' tp$setLikelihoodApproximator( approximatorVersion$AUTO_TUNING )
#'
#' # this is equivalent to:
#' tp$setLikelihoodApproximator( 0 )
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

#' tensorphylo debug modes
#'
#' @description
#' A list of debug modes for tensorphylo
#'
#' @details
#' Use these variables to specify debug modes
#' There are five options:
#' | **Setting** | **Index** | **Interpretation** |
#' |:--------|:-:|:--------------|
#' | `debugMode$DBG_NONE` _(default)_  | 0 | No debug info
#' | `debugMode$DBG_PRINT`             | 1 | Debug info printed to screen
#' | `debugMode$DBG_FILE`              | 2 | Debug info printed to file
#'
#' @examples
#' # create an empty tensorphylo object
#' tp <- new(TensorPhylo, 4)
#'
#' # specify a debug setting
#' tp$debugMode( debugMode$DBG_PRINT )
#'
#' # this is equivalent to:
#' tp$debugMode( 0 )
#' @name debugMode
#' @export
NULL

debugMode <- list(
  DBG_NONE=0,
  DBG_PRINT=1,
  DBG_FILE=2
)
