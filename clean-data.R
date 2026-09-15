source("functions.R")

# 1. Compute cartesian coordinates
# 2. Detrend and descale data
# 3. Plot data with trend (a chosen year is highlighted)
# 4. Plot data with scale (a chosen year is highlighted)

# Parse arguments
args <- OptionParser() |>
  add_option(c("-y", "--year_group"), type = "integer", default = 2025,
             help = "Year to highlight in time series plots") |>
  add_option(c("-w", "--half_width"), type = "integer", default = 30,
             help = "Window size = half_width *  2 for trend and scale") |>
  parse_args()

# For all the figures, the file is named according to the arguments
file_arg <- str_interp("_y-${year_group}_w-${half_width}", args)

#' Compute seasonal trend
#'
#' The trend function varies by day of the year 1,...,365 but is same for each
#' year.
#'
#' @param coord Matrix of wind vector coordinates. Each row represents one
#'   observations.
#' @param doy Vector where each element corresponds to the day of the year (DOY)
#'  of a coordinate. Values range from 1 to 365.
#' @param half_width Window size for each year is 2 * half_width.
#' @param n_iter Number of iterations for scale and scatter estimation.
#'
#' @returns A list giving trend for each DOY, scale for each DOY, estimated
#'   detrended and descaled observations, Mahalanobis distances,
#'   estimated generating variate and scatter for descaled and detrended
#'   observations.
fit_seasonal_window <- function(coord, doy, half_width, n_iter) {
  n_doy <- 365

  # Compute circular distances between DOY and DOY corresponding to an
  # observations
  gap <- abs(outer(seq_len(n_doy), doy, "-"))
  d_circular <- pmin(gap, n_doy - gap)

  # (i, j) element is true if jth observation lies in the window of day i
  w <- (d_circular <= half_width)
  day_window <- rowSums(w)
  if (any(day_window == 0)) {
    cli::cli_abort(c(
      "Widen the window",
      "x" = "Some day of the year has an empty window"
    ))
  }

  # Window means
  mu_doy <- (w %*% coord) / day_window
  mu_doy

  e <- coord - mu_doy[doy, ]
  s <- rep(1, nrow(coord)) # Initial guess for the scale
  for (i in seq_len(n_iter)) {
    sigma <- cov(e / s)
    sigma <- sigma / det(sigma)^(1 / ncol(coord))
    maha_dist <- sqrt(mahalanobis(e, center = FALSE, cov = sigma))
    s_doy <- (w %*% maha_dist) / day_window
    s <- s_doy[doy]
  }
  list(
    mu_doy = mu_doy,
    s_doy = s_doy,
    sigma = sigma,
    maha_dist = maha_dist,
    r = maha_dist / s,
    coord_0 = e / s
  )
}

# Delete station name column
# Rename columns for convenience
# Remove NA ("-" in the original data)
# Remove corrupted wind direction (there is one = 8348)
# Zero wind speeds have "a direction" but this does not matter in the
# transformation to cartesian coordinates
# Create a date column, get rid of hours (not needed)
wind <- readr::read_csv(
  "data/raw-wind.csv",
  col_types = "ciiitcc",
  col_names = TRUE
) |>
  select(-Havaintoasema) |>
  rename(
    year = Vuosi,
    month = Kuukausi,
    day = Päivä,
    time = `Aika [Paikallinen aika]`,
    speed = `Keskituulen nopeus [m/s]`,
    direction = `Tuulen suunnan keskiarvo [°]`
  ) |>
  mutate(
    speed = as.numeric(na_if(speed, "-")),
    direction = as.integer(na_if(direction, "-"))
  ) |>
  tidyr::drop_na() |>
  filter(direction >= 1  & direction <= 360) |>
  mutate(date = make_date(year, month, day)) |>
  select(date, speed, direction)
wind

# Transformation to cartesian coordinates
# Compute daily averages of the wind coordinates
# Include only days with at least 20 measurements
# For simplicity, delete leap days
wind_cartesian <- wind |>
  mutate(
    theta = direction * pi / 180,
    x = -speed * sin(theta),
    y = -speed * cos(theta)
  ) |>
  select(-(speed:theta)) |>
  summarise(
    x = mean(x),
    y = mean(y),
    n_hours = n(),
    .by = c(date)
  ) |>
  filter(n_hours >= 20) |>
  select(-n_hours) |>
  filter(!(month(date) == 2 & day(date) == 29))

coord <- wind_cartesian |>
  select(x, y) |>
  as.matrix()

# Compute day of the year (DOY) for each date. Days after the leap day are
# lagged so that DOY is in the range 1--365 for leap years also.
doy_leap <- 60
doy <- yday(wind_cartesian$date)
leap <- (year(wind_cartesian$date) %% 4 == 0)
doy <- ifelse(leap & doy >= doy_leap, doy - 1L, doy)

# Detrend and descale
estimates <- fit_seasonal_window(coord, doy, args$half_width, 5)

wind_cartesian <- wind_cartesian |>
  mutate(
    x_0 = estimates$coord_0[, 1],
    y_0 = estimates$coord_0[, 2],
    mu_x = estimates$mu_doy[doy, 1],
    mu_y = estimates$mu_doy[doy, 2],
    scale = estimates$s_doy[doy],
    r = estimates$r,
    maha_dist = estimates$maha_dist
  )
wind_cartesian

# Data for plotting
wind_plot <- wind_cartesian |>
  mutate(
    year_group = case_when(
      year(date) == args$year_group ~ "selected",
      TRUE ~ "other"
    )) |>
  mutate(md = as.Date(format(date, "2000-%m-%d"))) |>
  select(-date) |>
  tidyr::pivot_longer(
    -c(year_group, md),
    names_to = "label",
    values_to = "series"
  )

# Plot x-coordinate and trend
wind_plot |>
  filter(label == "x" |  label == "mu_x") |>
  ggplot(aes(md, series)) +
  geom_point(
    data = ~ filter(.x, label == "x" & year_group == "other"),
    colour = "grey40", alpha = 0.25, size = 1
  ) +
  geom_line(
    data = ~ filter(.x, label == "x" & year_group == "selected"),
    linetype = "longdash", linewidth = 0.8
  ) +
  geom_line(
    data = ~ filter(.x, label == "mu_x") |> distinct(md, series),
    linewidth = 1
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("x-coordinate") +
  xlab("Time") +
  theme_plot
ggsave(str_c("figures/x-series", file_arg, ".pdf"), dpi = 600)

# Plot y-coordinate and trend
wind_plot |>
  filter(label == "y" |  label == "mu_y") |>
  ggplot(aes(md, series)) +
  geom_point(
    data = ~ filter(.x, label == "y" & year_group == "other"),
    colour = "grey40", alpha = 0.25, size = 1
  ) +
  geom_line(
    data = ~ filter(.x, label == "y" & year_group == "selected"),
    linetype = "longdash", linewidth = 0.8
  ) +
  geom_line(
    data = ~ filter(.x, label == "mu_y") |> distinct(md, series),
    linewidth = 1
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("y-coordinate") +
  xlab("Time") +
  theme_plot
ggsave(str_c("figures/y-series", file_arg, ".pdf"), dpi = 600)

# Plot scale
wind_plot |>
  filter(label %in% c("maha_dist", "scale")) |>
  ggplot(aes(md, series)) +
  geom_point(
    data = ~ filter(.x, label == "maha_dist" & year_group == "other"),
    colour = "grey40", alpha = 0.25, size = 1
  ) +
  geom_line(
    data = ~ filter(.x, label == "maha_dist" & year_group == "selected"),
    linetype = "longdash", linewidth = 0.8
  ) +
  geom_line(
    data = ~ filter(.x, label == "scale") |> distinct(md, series),
    linewidth = 1
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  xlab("Time") +
  ylab("Scale") +
  theme_plot
ggsave(str_c("figures/scale", file_arg, ".pdf"), dpi = 600)

# Plot transformed data
wind_cartesian |>
  ggplot(aes(x = x_0, y = y_0)) +
  geom_point(colour = "grey40", alpha = 0.25, size = 1) +
  xlab("x-coordinate") +
  ylab("y-coordinate") +
  theme_plot
file_arg <- str_interp("_w-${half_width}", args)
ggsave(str_c("figures/scatter", file_arg, ".pdf"), dpi = 600)

# Write data
wind_cartesian |>
  readr::write_csv(str_c("data/cartesian-wind", file_arg, ".csv"))
