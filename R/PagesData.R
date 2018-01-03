## Copyright (C) 2015-2018 Thorsten Pohlert, GPL-3
#' @name PagesData
#' @docType data
#' @title Simulated data of Page (1955) as test-example
#' for change-point detection
#'
#' @description
#' Simulated data of Page (1955) as test-example for change-point
#' detection taken from Table 1 of Pettitt (1979)
#'
#' @usage data(PagesData)
#' 
#' @format a vector that contains 40 elements 
#'
#' @details
#'  According to the publication of Pettitt (1979), the series comprise a
#'  significant \eqn{p = 0.014} change-point at \eqn{i = 17}. The function
#'  \code{\link{pettitt.test}} computes the same U statistics as given by
#'  Pettitt (1979) in Table1, row 4.
#'
#' @references
#' Page, E. S. (1954), A test for a change in a parameter occuring at an
#' unknown point. \emph{Biometrika} 41, 100--114.
#'
#' Pettitt, A. N., (1979). A non-parametric approach to the change point
#' problem. \emph{Journal of the Royal Statistical Society Series C, Applied
#' Statistics} 28, 126--135.
#'
#' @seealso \code{\link{pettitt.test}}
#'
#' @examples
#' data(PagesData)
#' pettitt.test(PagesData)
#' @keywords datasets
NULL

