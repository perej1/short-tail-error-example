source("functions.R")

# 1. Compute quantile regions corresponding to windy and calm day
# 2. Plot region with days that are plusminus 30 days away from the
# windy/calm day

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
  select(date, x, y) |>
  mutate(
    maha_d = mahalanobis(
      x = bind_cols(x, y, .name_repair = ~ c("x", "y")),
      center = mu_max,
      cov = scale_max * sigma_est
    ),
    label = if_else(maha_d == max(maha_d), "max", "other")
  )


# Points close to calm day
day_min <- wind_min |> pull(date) |> yday()
d <- abs(yday(wind$date) - day_min)

wind_min_window <- wind |>
  mutate(in_window = (pmin(d, 365 - d) <= half_width)) |>
  filter(in_window == TRUE) |>
  select(date, x, y) |>
  mutate(
    maha_d = mahalanobis(
      x = bind_cols(x, y, .name_repair = ~ c("x", "y")),
      center = mu_min,
      cov = scale_min * sigma_est
    ),
    label = if_else(maha_d == max(maha_d), "max", "other")
  )


# Plotting
axislim <- c(-15, 15)
ggplot() +
  geom_path(data = region_max, aes(x = x, y = y, linetype = label)) +
  geom_point(
    data = wind_max_window,
    aes(x = x, y = y, shape = label, colour = label, size = label)
  ) +
  scale_linetype_manual(
    values = c("inner" = "solid", "outer" = "dashed")
  ) +
  scale_shape_manual(
    values = c("other" = 16, "max" = 17)
  ) +
  scale_size_manual(
    values = c("other" = 1, "max" = 4)
  ) +
  scale_colour_manual(
    values = c("other" = "#80808080", "max" = "black")
  ) +
  coord_equal(xlim = axislim, ylim = axislim) +
  xlab("x-coordinate") +
  ylab("y-coordinate") +
  theme_plot
file_arg <- str_interp("_w-${half_width}")
ggsave(str_c("figures/region_max", file_arg, ".pdf"), dpi = 600)

ggplot() +
  geom_path(data = region_min, aes(x = x, y = y, linetype = label)) +
  geom_point(
    data = wind_min_window,
    aes(x = x, y = y, shape = label, colour = label, size = label)
  ) +
  scale_linetype_manual(
    values = c("inner" = "solid", "outer" = "dashed")
  ) +
  scale_shape_manual(
    values = c("other" = 16, "max" = 17)
  ) +
  scale_size_manual(
    values = c("other" = 1, "max" = 4)
  ) +
  scale_colour_manual(
    values = c("other" = "#80808080", "max" = "black")
  ) +
  coord_equal(xlim = axislim, ylim = axislim) +
  xlab("x-coordinate") +
  ylab("y-coordinate") +
  theme_plot
ggsave(str_c("figures/region_min", file_arg, ".pdf"), dpi = 600)

# Print extreme obs wrt windy/calm day
show_windy_date <- wind_max |>
  pull(date) |>
  format("xxxx-%m-%d")
show_windy_max_date <- wind_max_window |>
  filter(label == "max") |>
  pull(date)
cli::cli_alert_info("Extreme observation on the date {show_windy_date}: {show_windy_max_date}")

show_calm_date <- wind_min |>
  pull(date) |>
  format("xxxx-%m-%d")
show_calm_min_date <- wind_min_window |>
  filter(label == "max") |>
  pull(date)
cli::cli_alert_info("Extreme observation on the date {show_calm_date}: {show_calm_min_date}")
