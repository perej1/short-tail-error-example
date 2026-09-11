library(optparse)
library(stringr)

# Parse arguments
args <- OptionParser() |>
  add_option(c("-w", "--half_width"), type = "integer", default = 30,
             help = "Window size = half_width *  2 for trend and scale") |>
  parse_args()

args

wind <- readr::read_csv(
  str_interp("data/cartesian-wind_w-${half_width}.csv", args),
  col_names = TRUE,
  col_types = str_flatten(c("D", rep("d", 9)))
)
