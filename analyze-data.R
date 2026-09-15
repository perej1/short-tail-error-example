source("functions.R")


# Parse arguments
args <- OptionParser() |>
  add_option(c("-w", "--half_width"), type = "integer", default = 30,
             help = "Window size = half_width *  2 for trend and scale") |>
  parse_args()

wind <- readr::read_csv(
  str_interp("data/cartesian-wind_w-${half_width}.csv", args),
  col_names = TRUE,
  col_types = str_flatten(c("D", rep("d", 9)))
)

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

# TODO: Estimate regions and extreme value index
