# short-tail-error-example
Illustrations for a paper about inference with contaminated (short-tailed) extremes

## Original data

The original data set `data/raw-wind.csv` is sourced from the [Finnish
Meteorological Institute](https://en.ilmatieteenlaitos.fi/open-data), and it is
provided under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
license:

- Retrieved on 3 September 2026

- Data range is 20 March 2003--31 December 2025

- Variables are hourly average windspeed and hourly average wind direction. Wind
  speed is measured in meters per second. Direction is measured in degrees
  (0/360 &rarr; wind from north to south, 90 &rarr; wind from east to west, 180
  &rarr; wind from south to north, 270 &rarr; wind from west to east).

- The bivariate wind time series is measured in the *Helsinki-Vantaan
  Lentoasema* station.

## Cleaned data

We translate the wind coordinates from polar to cartesian. That is, the
transformed coordinates are wind vectors, where the length of the vector
represents speed and directions are given on the below table.

(x,y)-direction  | Wind direction |
| ------------- | ------------- |
| Positive y-direction  | Wind from south to north |
| Negative y-direction  | Wind from north to south |
| Positive x-direction  | Wind from west to east |
| Negative x-direction  | Wind from east to west |

After transformation to cartesian coordinates, the hourly data is averaged to
daily level. The cleaned data is in the file `cartesian-wind.csv`.