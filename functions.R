library(lubridate)
library(ggplot2)
library(optparse)
library(stringr)
library(dplyr)

# Global constants and functions

# My plot theme for ggplot2
theme_plot <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    text       = element_text(size = 16),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),

    panel.grid.major = element_line(color = scales::alpha("black", 0.2)),
    panel.grid.minor = element_line(color = scales::alpha("black", 0.1)),

    legend.position = "none"
  )


#' Convert a spherical coordinate to a Cartesian coordinate
#'
#' @param radius Double, radius of the point.
#' @param theta Double vector, d - 1 angular coordinates.
#'
#' @return Double vector, d coordinates giving the location of the point in
#'   space.
spherical_to_cartesian <- function(radius, theta) {
  d <- length(theta) + 1
  cartesian <- rep(NA, d)
  cartesian[1] <- cos(theta[1])
  cartesian[d] <- prod(sin(theta))

  if (d > 2) {
    for (i in 2:(d - 1)) {
      cartesian[i] <- prod(sin(theta[1:(i - 1)])) * cos(theta[i])
    }
  }
  radius * cartesian
}


#' Generate m points from a (d - 1)-sphere
#'
#' @param d Integer, dimension of the ball, must be equal to or greater than 2.
#' @param m_angle Integer, number of points to return.
#'
#' @return List of two double matrices, the first gives the ball in spherical
#'   and the second in Cartesian coordinates. One row represents one point.
get_ball_mesh <- function(d, m_angle) {
  m <- ceiling(m_angle^(1 / (d - 1)))
  angle_list <- vector("list", d - 1)
  angle_list[[d - 1]] <- seq(0, 2 * pi, length.out = m)
  if (d > 2) {
    for (i in 1:(d - 2)) {
      angle_list[[i]] <- seq(0, pi, length.out = m)
    }
  }
  spherical <- as.matrix(expand.grid(angle_list))
  list(spherical = spherical,
       cartesian = t(apply(spherical, 1, spherical_to_cartesian, r = 1)),
       m_effective = nrow(spherical))
}


#' Compute square root of a positive definite matrix
#'
#' More precisely, function computes matrix lambda
#' s.t. sigma = lambda %*% lambda.
#'
#' @param sigma Double matrix, positive definite matrix.
#'
#' @return Double matrix, square root of the matrix.
sqrtmat <- function(sigma) {
  eigenval <- eigen(sigma)$values
  if (any(eigenval <= 0) || any(sigma != t(sigma))) {
    rlang::abort("`sigma` must be a symmetric positive definite matrix.")
  }
  eigenvec <- eigen(sigma)$vectors
  eigenvec %*% diag(eigenval^0.5) %*% t(eigenvec)
}


#' Helper function for moment-based estimation
#'
#' @param data_inc Double data vector in increasing order.
#' @param n Integer, sample size.
#' @param k Integer, threshold for tail observations.
#' @param l Integer, exponent.
#'
#' @returns Double.
moment_component <- function(data_inc, n, k, l) {
  mean(((log(data_inc[(n - k):n]) - log(data_inc[n - k]))[-1])^l)
}


#' Estimate negative part of the extreme value index
#'
#' @param data_inc Double data vector in increasing order.
#' @param n Integer, sample size.
#' @param k Integer, threshold for tail observations.
#'
#' @returns Double, the estimate.
estimate_gamma_minus <- function(data_inc, n, k) {
  up <- moment_component(data_inc, n, k, 1)^2
  low <- moment_component(data_inc, n, k, 2)
  1 - 0.5 * (1 - up / low)^(-1)
}


#' Estimate the extreme value index
#'
#' @param data_inc Double data vector in increasing order.
#' @param n Integer, sample size.
#' @param k Integer, threshold for tail observations.
#'
#' @returns Double, the estimate.
estimate_gamma <- function(data_inc, n, k) {
  moment_component(data_inc, n, k, 1) + estimate_gamma_minus(data_inc, n, k)
}


#' Estimate (1-p)-quantile
#'
#' @param data_inc Double data vector in increasing order.
#' @param n Integer, sample size.
#' @param k Integer, threshold for tail observations.
#' @param p Double, probability corresponding to the (1-p)-quantile
#'
#' @returns Double, the estimate.
estimate_quantile <- function(data_inc, n, k, p) {
  sigma_m <- data_inc[n - k] * moment_component(data_inc, n, k, 1) *
    (1 - estimate_gamma_minus(data_inc, n, k))
  up <- (k / (n * p))^estimate_gamma(data_inc, n, k) - 1
  data_inc[n - k] + sigma_m * up / estimate_gamma(data_inc, n, k)
}


#' Estimate elliptical extreme quantile region
#'
#' @param data Double matrix, each row corresponds to one observation.
#' @param mu_est Double vector, location estimate.
#' @param sigma_est Double matrix, estimate of the shape.
#' @param k Integer, threshold for tail observations.
#' @param p Double, probability corresponding to the (1-p)-quantile.
#' @param m_angle Integer, number of points to return from the boundary of the
#'   estimate.
#'
#' @returns List of the following objects: 1. A matrix giving m_angle points
#'   from the boundary of the quantile region, 2. estimated (1 - p)-quantile of
#'   the generating variate, 3. estimated extreme value index for the elliptical
#'   distribution.
estimate_qregion <- function(data, mu_est, sigma_est, k, p, m_angle) {
  n <- nrow(data)
  d <- ncol(data)
  w <- get_ball_mesh(d, m_angle)$cartesian

  # Center data
  data <- sweep(data, 2, mu_est, "-")

  # Approximate generating variate
  radius <- sqrt(stats::mahalanobis(data, FALSE, sigma_est, inverted = FALSE))
  radius_inc <- sort(radius, decreasing = FALSE)

  # Estimate extreme value index
  gamma_est <- estimate_gamma(radius_inc, n, k)

  # Estimate extreme quantile of the generating variate
  r_hat <- estimate_quantile(radius_inc, n, k, p)

  # Estimate extreme quantile region
  lambda <- sqrtmat(r_hat^2 * sigma_est)
  region <- sweep(w %*% t(lambda), 2, mu_est, "+")
  list(region = region, r_hat = r_hat, gamma_est = gamma_est)
}
