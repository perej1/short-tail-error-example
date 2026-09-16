source("functions.R")

# Global constants
half_width <- 30
k <- 1500
m_angle <- 1000

wind <- readr::read_csv(
  str_interp("data/cartesian-wind_w-${half_width}.csv"),
  col_names = TRUE,
  col_types = str_flatten(c("D", rep("d", 9)))
)


# Estimate quantile region for standardized data
data <- wind |>
  select(x_0, y_0) |>
  as.matrix()
mu_est <- c(0, 0)
sigma_est <- cov(data)
p <- c(0.01, 1 / (2 * nrow(data)))


inner <- estimate_qregion(data, mu_est, sigma_est, k, p[1], m_angle)
outer <- estimate_qregion(data, mu_est, sigma_est, k, p[2], m_angle)
region_inner <- tibble::as_tibble(inner$region, .name_repair = ~ c("x", "y"))
region_outer <- tibble::as_tibble(outer$region, .name_repair = ~ c("x", "y"))

region <- bind_rows(region_inner, region_outer) |>
  mutate(label = c(rep("inner", m_angle), rep("outer", m_angle)))

# Choose windy day
wind_max <- wind |>
  slice_max(scale, n = 1) |>
  slice_max(r, n = 1)

scale_max <- wind_max |> pull(scale)
mu_max <- wind_max |>
  select(mu_x, mu_y) |>
  as.vector() |>
  unlist(use.names = FALSE)

region_max <- region |>
  mutate(
    x = scale_max * x + mu_max[1],
    y = scale_max * y + mu_max[2]
  )

# Choose calm day
wind_min <- wind |>
  slice_min(scale, n = 1) |>
  slice_min(r, n = 1)
scale_min <- wind_min |> pull(scale)
mu_min <- wind_min |>
  select(mu_x, mu_y) |>
  as.vector() |>
  unlist(use.names = FALSE)

region_min <- region |>
  mutate(
    x = scale_min * x + mu_min[1],
    y = scale_min * y + mu_min[2]
  )


# Points close to windy day
day_max <- wind_max |> pull(date) |> yday()
d <- abs(yday(wind$date) - day_max)

wind_max_window <- wind |>
  mutate(in_window = (pmin(d, 365 - d) <= half_width)) |>
  filter(in_window == TRUE) |>
  select(date, x, y)

# Points close to calm day
day_min <- wind_min |> pull(date) |> yday()
d <- abs(yday(wind$date) - day_min)

wind_min_window <- wind |>
  mutate(in_window = (pmin(d, 365 - d) <= half_width)) |>
  filter(in_window == TRUE) |>
  select(date, x, y)

# Plotting
axislim <- c(-15, 15)
ggplot() +
  geom_path(data = region_max, aes(x = x, y = y, linetype = label)) +
  geom_point(
    data = wind_max_window,
    aes(x = x, y = y),
    colour = "grey40",
    alpha = 0.5,
    size = 1.5
  ) +
  scale_linetype_manual(
    values = c("inner" = "solid", "outer" = "dashed")
  ) +
  coord_equal(xlim = axislim, ylim = axislim) +
  theme_plot

ggplot() +
  geom_path(data = region_min, aes(x = x, y = y, linetype = label)) +
  geom_point(
    data = wind_min_window,
    aes(x = x, y = y),
    colour = "grey40",
    alpha = 0.5,
    size = 1.5
  ) +
  scale_linetype_manual(
    values = c("inner" = "solid", "outer" = "dashed")
  ) +
  coord_equal(xlim = axislim, ylim = axislim) +
  theme_plot
