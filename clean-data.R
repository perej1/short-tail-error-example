library(dplyr)
library(lubridate)

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
) %>%
  select(-Havaintoasema) %>%
  rename(
    year = Vuosi,
    month = Kuukausi,
    day = Päivä,
    time = `Aika [Paikallinen aika]`,
    speed = `Keskituulen nopeus [m/s]`,
    direction = `Tuulen suunnan keskiarvo [°]`
) %>%
  mutate(
    speed = as.numeric(na_if(speed, "-")),
    direction = as.integer(na_if(direction, "-"))
) %>%
  tidyr::drop_na() %>%
  filter(direction >= 1  & direction <= 360) %>%
  mutate(date = make_date(year, month, day)) %>%
  select(date, speed, direction)
wind

# Transformation to cartesian coordinates
# Compute daily averages of the wind coordinates
# Include only days with at least 20 measurements
wind_cartesian <- wind %>%
  mutate(
    theta = direction * pi / 180,
    x = -speed * sin(theta),
    y = -speed * cos(theta)
) %>%
  select(-(speed:theta)) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    n_hours = n(),
    .by = c(date)
  ) %>%
  filter(n_hours >= 20) %>%
  select(-n_hours)
wind_cartesian

# Write data
wind_cartesian %>%
  readr::write_csv("data/cartesian-wind.csv")
