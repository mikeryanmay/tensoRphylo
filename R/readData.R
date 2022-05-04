#' Read Nexus formatted data
#'
#' Read data in the Nexus file format for tensorphylo.
#' In principle, works with any nexus datatype and symbols (I haven't done extensive testing), and handles ambiguity/missing data.
#'
#' @param file path to a nexus file with one DATA block, that contains a single discrete character.
#' @return A matrix formatted for input to tensorphylo. Named rows correspond to species, columns correspond to character states. A value of 1 indicates that the character is consistent with the observed data.
#'
#' @name readNexusData
#' @export
NULL

#' Read character data in a delimited text file for tensorphylo.
#'
#' @details
#' Read data in a delimited text file.
#' Useful for cases where nexus doesn't accommodate your dataset (e.g., discrete characters with lots of states, like chromosome numbers).
#' The file may be comma- or tab-delimited, with the first column corresponding to the species name, and the second column corresponding to the character state(s).
#' States should be encoded as numbers between 0 and `nstates - 1` (i.e., the first state is 0).
#' Ambiguous data can be encoded by listing all consistent states, separated by `/`. For example, `2/3/7` indicates the species is possibly in state 2, 3 or 7.
#' Question marks, `?`, or dashes ,`-`, indicate missing data.
#'
#' @param file path to a delimited text file (see details for format).
#' @param delim the delimiter, either ",", or "/".
#' @param nstates the number of states.
#' @return A matrix formatted for input to tensorphylo. Named rows correspond to species, columns correspond to character states. A value of 1 indicates that the character is consistent with the observed data.
#'
#' @name readDelimitedData
#' @export
NULL
