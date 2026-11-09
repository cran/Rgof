#' power_study
#' 
#' the results of the included power studies
#' @format ## 'power_study'
#' A list with parts cont and disc, each a list of matrices with powers (case studies
#' in the rows, methods in the columns). attr(power_study, "Null") has the type I
#' errors, a list with matrices cont and disc.
"power_study"

#' power_study_hybrid
#' 
#' the results of the included power studies for the hybrid tests (hybrid_power),
#' with a Monte Carlo data set of the same size as the data (nMC=1)
#' @format ## 'power_study_hybrid'
#' As power_study, including attr(power_study_hybrid, "Null").
"power_study_hybrid"

#' power_study_hybrid_nMC5
#' 
#' the results of the included power studies for the hybrid tests (hybrid_power),
#' with a Monte Carlo data set five times the size of the data (nMC=5)
#' @format ## 'power_study_hybrid_nMC5'
#' As power_study, including attr(power_study_hybrid_nMC5, "Null").
"power_study_hybrid_nMC5"

#' pvaluecdf
#' 
#' the info needed to draw a graph
#' @format ## 'pvaluecdf'
#' A matrix
"pvaluecdf"

#' case_studies_sample_sizes
#' 
#' sample sizes for the various case studies
#' @format ## 'case_studies_sample_sizes'
#' A list
"case_studies_sample_sizes"

#' funs_list
#' 
#' function definitions for case studies
#' @format ## 'funs_list'
#' A list of functions
"funs_list"

#' power.uniform.linear
#' 
#' matrix with power values, for vignette
#' @format ## 'power.uniform.linear'
#' matrix with power values, for vignette
"power.uniform.linear"
