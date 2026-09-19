source("functions.R")

# Boxplots of normalized errors in various scenarios

#' Generate sample from an elliptical distribution
#'
#' Location is zero, scatter is identity and generating variate is beta
#' distributed.
#'
#' @param n Integer, sample size.
#' @param gamma Extreme value index for the generating variate.
#'
#' @returns List of two elements: i) tibble, where each row is one observation
#' ii) real sample from the generating variate.
sample_elliptical <- function(n, gamma) {
  rsample <- rbeta(n, 1, 1 / abs(gamma))
  obs <- tibble(
    x = rnorm(n),
    y = rnorm(n)
  ) |>
    mutate(
      r = sqrt(x^2 + y^2),
      x = x / r * rsample,
      y = y / r * rsample
    ) |>
    select(-r)
  list(obs = obs, rsample = rsample)
}


#' Normalization for the errors
#'
#' @param n Integer, sample size.
#' @param k Integer, threshold for tail observations.
#' @param p Double, probability corresponding to the (1-p)-quantile.
#' @param gamma Integer, extreme value index.
#'
#' @returns Double giving the divisor for the errors.
standardize <- function(n, k, p, gamma) {
  x1 <- n / k
  x2 <- k / (n * p)
  q1 <- (x2^gamma * log(x2)) / gamma
  q2 <- (x2^gamma - 1) / gamma^2
  q <- q1 - q2
  a <- abs(gamma) * x1^gamma
  q * a
}


#' Compute estimation errors for various quantile estimators
#'
#' @param m Number of repetitions.
#' @param seed Integer, random seed for sample generation.
#' @param n Integer, sample size.
#' @param gamma Extreme value index for the generating variate.
#' @param p Double, probability corresponding to the (1-p)-quantile.
#' @param k Integer, threshold for tail observations.
#' @param j Integer, simulation round.
#' @param n_rounds Integer, number of simulation rounds.
#'
#' @returns Tibble, each column has m errors corresponding to an estimator.
simulate <- function(m, seed, n, gamma, p, k, j, n_rounds) {
  cli::cli_alert_info("Round {j} / {n_rounds}")
  set.seed(seed)
  error_er <- rep(NA, m)
  error_ee <- rep(NA, m)
  error_nr <- rep(NA, m)
  error_ne <- rep(NA, m)
  for (i in 1:m) {
    # Generate sample
    elliptical <- sample_elliptical(n, gamma)
    obs <- as.matrix(elliptical$obs)
    rsample_inc <- sort(elliptical$rsample, decreasing = FALSE)

    # Approximations of the generating variate
    mu <- colMeans(obs)
    sigma <- cov(obs)
    sigma <- sigma / det(sigma)^(1 / ncol(obs))
    radius <- sqrt(stats::mahalanobis(obs, mu, sigma))
    radius_inc <- sort(radius, decreasing = FALSE)

    # Standardized errors
    q_real <- qbeta(1 - p, 1, 1 / abs(gamma))
    error_er[i] <- (estimate_quantile(rsample_inc, n, k, p) - q_real) / standardize(n, k, p, gamma)
    error_ee[i] <- (estimate_quantile(radius_inc, n, k, p) - q_real) / standardize(n, k, p, gamma)
    error_nr[i] <- (rsample_inc[n] - q_real) / standardize(n, k, p, gamma)
    error_ne[i] <- (radius_inc[n] - q_real) / standardize(n, k, p, gamma)
  }
  tibble(
    error_er = error_er,
    error_ee = error_ee,
    error_nr = error_nr,
    error_ne = error_ne
  )
}

# parameters
params <- tidyr::expand_grid(
  m = 500,
  seed = 123,
  n = 10^(4:6),
  gamma = c(-0.25, -0.75),
  k_exp = c(0.3, 0.6)
) |>
  mutate(
    p = 1 / n,
    k = floor(n^k_exp)
  )

n_rounds <- nrow(params)
params <- params |>
  mutate(
    j = 1:n_rounds,
    n_rounds = n_rounds
  )

# Perform simulations
results <- params |>
  mutate(
    simulation = purrr::pmap(
      list(m, seed, n, gamma, p, k, j, n_rounds),
      simulate
    )
  )


# Plotting
results |>
  filter(gamma == -0.25 & k_exp == 0.3) |>
  tidyr::unnest(simulation) |>
  tidyr::pivot_longer(
    cols = c(error_er, error_ee, error_nr, error_ne),
    names_to = "method",
    values_to = "error"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("error_er", "error_ee", "error_nr", "error_ne")
    )
  ) |>
  ggplot(aes(x = factor(n), y = error, fill = method)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      error_er = "steelblue",
      error_ee = "tomato",
      error_nr = "goldenrod",
      error_ne = "white"
    )
  ) +
  geom_hline(yintercept = 0) +
  xlab("Sample size") +
  ylab("Normalized error") +
  theme_plot
ggsave("figures/small_gamma_small_k.pdf", dpi = 600)

results |>
  filter(gamma == -0.25 & k_exp == 0.6) |>
  tidyr::unnest(simulation) |>
  tidyr::pivot_longer(
    cols = c(error_er, error_ee, error_nr, error_ne),
    names_to = "method",
    values_to = "error"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("error_er", "error_ee", "error_nr", "error_ne")
    )
  ) |>
  ggplot(aes(x = factor(n), y = error, fill = method)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      error_er = "steelblue",
      error_ee = "tomato",
      error_nr = "goldenrod",
      error_ne = "white"
    )
  ) +
  geom_hline(yintercept = 0) +
  xlab("Sample size") +
  ylab("Normalized error") +
  theme_plot
ggsave("figures/small_gamma_large_k.pdf", dpi = 600)

results |>
  filter(gamma == -0.75 & k_exp == 0.3) |>
  tidyr::unnest(simulation) |>
  tidyr::pivot_longer(
    cols = c(error_er, error_ee, error_nr, error_ne),
    names_to = "method",
    values_to = "error"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("error_er", "error_ee", "error_nr", "error_ne")
    )
  ) |>
  ggplot(aes(x = factor(n), y = error, fill = method)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      error_er = "steelblue",
      error_ee = "tomato",
      error_nr = "goldenrod",
      error_ne = "white"
    )
  ) +
  geom_hline(yintercept = 0) +
  xlab("Sample size") +
  ylab("Normalized error") +
  theme_plot
ggsave("figures/large_gamma_small_k.pdf", dpi = 600)

results |>
  filter(gamma == -0.75 & k_exp == 0.6) |>
  tidyr::unnest(simulation) |>
  tidyr::pivot_longer(
    cols = c(error_er, error_ee, error_nr, error_ne),
    names_to = "method",
    values_to = "error"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("error_er", "error_ee", "error_nr", "error_ne")
    )
  ) |>
  ggplot(aes(x = factor(n), y = error, fill = method)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      error_er = "steelblue",
      error_ee = "tomato",
      error_nr = "goldenrod",
      error_ne = "white"
    )
  ) +
  geom_hline(yintercept = 0) +
  xlab("Sample size") +
  ylab("Normalized error") +
  theme_plot
ggsave("figures/large_gamma_large_k.pdf", dpi = 600)
