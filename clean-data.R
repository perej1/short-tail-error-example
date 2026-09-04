library(dplyr)
library(lubridate)
library(ggplot2)
source("theme-plot.R")


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
#'
#' @returns Vector giving the trend for each DOY
fit_seasonal_window <- function(coord, doy, half_width) {
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
wind_cartesian

coord <- wind_cartesian |>
  select(x, y) |>
  as.matrix()

# Compute day of the year (DOY) for each date. Days after the leap day are
# lagged so that DOY is in the range 1--365 for leap years also.
doy_leap <- 60
doy <- yday(wind_cartesian$date)
leap <- (year(wind_cartesian$date) %% 4 == 0)
doy <- ifelse(leap & doy >= doy_leap, doy - 1L, doy)

# Compute trend
half_width <- 30
mu_doy <- fit_seasonal_window(coord, doy, half_width)

wind_cartesian <- wind_cartesian |>
  mutate(mu_x = mu_doy[doy, 1], mu_y = mu_doy[doy, 2])

# Time series plots of coordinates for a chosen period
wind_plot <- wind_cartesian |>
  filter(year(date) %in% 2024:2025) |>
  tidyr::pivot_longer(!date, names_to = "label", values_to = "series")

wind_plot |>
  filter(label %in% c("x", "mu_x")) |>
  ggplot(aes(x = date, y = series, linetype = label)) +
  geom_line() +
  xlab("Date") +
  ylab("x-coordinate") +
  scale_linetype_manual(
    values = c("x" = "dotted",  "mu_x" = "solid")
  )  +
  theme_plot
ggsave("figures/x-series.pdf", dpi = 600)

wind_plot |>
  filter(label %in% c("y", "mu_y")) |>
  ggplot(aes(x = date, y = series, linetype = label)) +
  geom_line() +
  xlab("Date") +
  ylab("y-coordinate") +
  scale_linetype_manual(
    values = c("y" = "dotted",  "mu_y" = "solid")
  )  +
  theme_plot
ggsave("figures/y-series.pdf", dpi = 600)

# Write data
wind_cartesian |>
  readr::write_csv("data/cartesian-wind.csv")
