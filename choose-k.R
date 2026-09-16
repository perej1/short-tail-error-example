source("functions.R")

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

# Compute approximations of the generating variate
data <- wind |>
  select(x_0, y_0) |>
  as.matrix()

radius <- sqrt(stats::mahalanobis(data, FALSE, cov(data), inverted = FALSE))
radius_inc <- sort(radius, decreasing = FALSE)

gamma_for_k <- tibble::tibble(
  k = 100:3000,
  gamma = purrr::map_dbl(k, ~ estimate_gamma(radius_inc, nrow(data), .))
)
readr::write_csv(
  gamma_for_k,
  str_c("results/k-cartesian-wind", file_arg, ".csv")
)

ggplot(gamma_for_k, aes(x = k, y = gamma)) +
  geom_line() +
  ylab("Estimate of the extreme value index") +
  theme_plot
ggsave(str_c("figures/k", file_arg, ".pdf"), dpi = 600)

gamma_neg <- all(gamma_for_k$gamma < 0)
gamma_neg
cli::cli_alert_info("All estimates of gamma negative? {gamma_neg}")
