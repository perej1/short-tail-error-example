source("functions.R")


# Parse arguments
args <- OptionParser() |>
  add_option(c("-w", "--half_width"), type = "integer", default = 30,
             help = "Window size = half_width *  2 for trend and scale") |>
  add_option(c("-k", "--k"), type = "integer", default = 10,
             help = "Threshold for the number of tail observations") |>
  parse_args()

wind <- readr::read_csv(
  str_interp("data/cartesian-wind_w-${half_width}.csv", args),
  col_names = TRUE,
  col_types = str_flatten(c("D", rep("d", 9)))
)

# TODO: Estimate regions and extreme value index

# mu <- as.vector(unlist(wind[ind, c(6, 7)]))
# mu
# 
# scale <- as.vector(unlist(wind[ind, c(8)]))
# scale
