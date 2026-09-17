source("functions.R")

# 1. Compute autocorrelation for original and thinned data
# 2. Compute gamma estimates for various values of k for original and thinned
# data

# Parse arguments
args <- OptionParser() |>
  add_option(c("-w", "--half_width"), type = "integer", default = 15,
             help = "Window size = half_width *  2 for trend and scale") |>
  parse_args()

wind <- readr::read_csv(
  str_interp("data/cartesian-wind_w-${half_width}.csv", args),
  col_names = TRUE,
  col_types = str_flatten(c("D", rep("d", 9)))
)

file_arg <- str_interp("_w-${half_width}", args)

# Thin data
wind_thin <- wind |>
  slice(seq(1, n(), by = 7))

# Compute autocorrelations
tibble::tibble(
  x_0 = as.vector(acf(wind$x_0, plot = FALSE, lag.max = 30)$acf),
  y_0 = as.vector(acf(wind$y_0, plot = FALSE, lag.max = 30)$acf),
  x_0_thin = as.vector(acf(wind_thin$x_0, plot = FALSE, lag.max = 30)$acf),
  y_0_thin = as.vector(acf(wind_thin$y_0, plot = FALSE, lag.max = 30)$acf)
) |>
  readr::write_csv(str_interp("results/acf_w-${half_width}.csv", args))

# Compute approximations of the generating variate
data <- wind |>
  select(x_0, y_0) |>
  as.matrix()

data_thin <- wind_thin |>
  select(x_0, y_0) |>
  as.matrix()

radius <- sqrt(stats::mahalanobis(data, FALSE, cov(data), inverted = FALSE))
radius_inc <- sort(radius, decreasing = FALSE)

radius_thin <- sqrt(stats::mahalanobis(data_thin, FALSE, cov(data_thin),
                                       inverted = FALSE))
radius_inc_thin <- sort(radius_thin, decreasing = FALSE)

gamma_for_k <- tibble::tibble(
  k = 100:3000,
  gamma = purrr::map_dbl(k, ~ estimate_gamma(radius_inc, nrow(data), .))
)
readr::write_csv(
  gamma_for_k,
  str_c("results/k-cartesian-wind", file_arg, ".csv")
)

gamma_for_k_thin <- tibble::tibble(
  k = 50:400,
  gamma = purrr::map_dbl(k, ~ estimate_gamma(radius_inc_thin, nrow(data_thin), .))
)
readr::write_csv(
  gamma_for_k,
  str_c("results/k-cartesian-wind-thin", file_arg, ".csv")
)


ggplot(gamma_for_k, aes(x = k, y = gamma)) +
  geom_line() +
  ylab("Estimate of the extreme value index") +
  theme_plot
ggsave(str_c("figures/k", file_arg, ".pdf"), dpi = 600)

gmin <- round(min(gamma_for_k$gamma), 2)
gmax <- round(max(gamma_for_k$gamma), 2)
cli::cli_alert_info("Range of gamma estimates ({gmin}, {gmax})")


ggplot(gamma_for_k_thin, aes(x = k, y = gamma)) +
  geom_line() +
  ylab("Estimate of the extreme value index") +
  theme_plot
ggsave(str_c("figures/k_thin", file_arg, ".pdf"), dpi = 600)

gmin <- round(min(gamma_for_k_thin$gamma), 2)
gmax <- round(max(gamma_for_k_thin$gamma), 2)
cli::cli_alert_info("Range of gamma estimates (thinned): ({gmin}, {gmax})")
