#' Graphical Horseshoe / Horseshoe-Lasso (RT) Sampler (C++ backend)
#'
#' @description
#' High-performance Reverse-Telescoping (RT) MCMC sampler for sparse precision
#' matrices under the Graphical Horseshoe (\code{"GHS"}) or Graphical Horseshoe Lasso
#' (\code{"GHSL"}) priors. This is a thin R wrapper to a C++ implementation exported
#' via Rcpp.
#'
#' @param Y Numeric matrix \eqn{n \times p}. Each column is a variable.
#' @param M Integer. Total MCMC iterations.
#' @param burnin Integer. Burn-in iterations to discard (\eqn{< M}).
#' @param seed Integer. RNG seed used for both \code{std::mt19937} and Armadillo RNG.
#' @param prior Character. One of \code{"GHS"} or \code{"GHSL"}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{Omega_save}: \eqn{p \times p \times (M - burnin)} posterior samples of the precision matrix.
#'   \item \code{Lambda2_save}: \eqn{p \times p \times (M - burnin)} local shrinkage \eqn{\lambda^2}.
#'   \item \code{tau2_save}: length-\eqn{(M - burnin)} vector of global parameter \eqn{\tau^2}.
#' }
#'
#' @export
#' @useDynLib RTsampler, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom GIGrvg rgig
#' @importFrom Matrix nearPD
#'
#' @name RT_sampler_cpp
#' @rdname RT_sampler_cpp
NULL
