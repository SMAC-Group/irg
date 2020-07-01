#' Raw measurements of two gene expressions of Arabidopsis thaliana
#'
#' A dataset containing three repeated measurements of two gene expressions (root and shoot respectively) from the fine scale time-series
#' RNA-Seq of shoot and root responses to Nitrogen supply in Arabidopsis thaliana from the Gene Expression Omnibus.
#'
#' @format A matrix with 30 rows and 3 columns (variables):
#' \describe{
#'    \item{root}{expression in the root}
#'    \item{shoot}{expression in the shoot}
#'    \item{time}{time from the beginning of the experiment (in minutes) at which measurements were taken}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97500}
"raw_data"

#' Detrended and standardized measurements of two gene expressions of Arabidopsis thaliana
#'
#' A list containing measurements of two pre-processed gene expressions (root and shoot respectively) from the fine scale time-series
#' RNA-Seq of shoot and root responses to Nitrogen supply in Arabidopsis thaliana from the Gene Expression Omnibus (whose raw measurements can be found in \code{raw_data}).
#'
#' @format A list with the following elements:
#' \describe{
#'    \item{root}{pre-processed expression in the root}
#'    \item{shoot}{pre-processed expression in the shoot}
#' }
#'
"signals"

#' Pre-processed measurements of gene expressions of Arabidopsis thaliana
#'
#' A list containing datasets of gene expressions measured over time in the root and shoot of Arabidopsis thaliana respectively. These were pre-processed
#' from the fine scale time-series RNA-Seq of shoot and root responses to Nitrogen supply in Arabidopsis thaliana from the Gene Expression Omnibus.
#'
#' @format A list with 2 elements:
#' \describe{
#'    \item{roots}{Data frame containing pre-processed gene expression measurements over time in the root (rows are different expressions and columns are measurement times)}
#'    \item{shoot}{Data frame containing pre-processed gene expression measurements over time in the shoot (rows are different expressions and columns are measurement times)}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97500}
"arabidopsis"

#' @details The data sets in this package may change at a moments notice.
"_PACKAGE"
