# short-tail-error-example

1. Empirical example about wind extremes. More precisely, elliptical extreme
quantile regions are computed for wind vectors.

2. Small simulation study that is independent of the empirical example.

## Empirical example

### Original data

The original data set `data/raw-wind.csv` is sourced from the [Finnish
Meteorological Institute](https://www.ilmatieteenlaitos.fi/havaintojen-lataus),
and it is provided under the [CC BY
4.0](https://creativecommons.org/licenses/by/4.0/) license:

- Retrieved on 3 September 2026

- Data range is 20 March 2003--31 December 2025

- Variables are hourly average windspeed and hourly average wind direction. Wind
  speed is measured in meters per second. Direction is measured in degrees
  (0/360 &rarr; wind from north to south, 90 &rarr; wind from east to west, 180
  &rarr; wind from south to north, 270 &rarr; wind from west to east).

- The bivariate wind time series is measured in the *Helsinki-Vantaan
  Lentoasema* station.

### Cleaned data

All the cleaning steps are done by the script `clean-data.R` First, we translate
the wind coordinates from polar to cartesian. That is, the transformed
coordinates are wind vectors, where the length of the vector represents speed
and directions are given by the below table.

(x,y)-direction  | Wind direction |
| ------------- | ------------- |
| Positive y-direction  | Wind from south to north |
| Negative y-direction  | Wind from north to south |
| Positive x-direction  | Wind from west to east |
| Negative x-direction  | Wind from east to west |

After the transformation to cartesian coordinates, the hourly data is averaged
to daily level for both coordinates (columns `x`  and `y`) on files of the
format `data/cartesian-wind_w-<half_width>.csv`.

For the Cartesian coordinates we suppose the model $(x_i,y_i)^T = \mu(t) + s(t)
\cdot (x_i^{(0)}, y_i^{(0)})^T$, where

- $\mu:\mathbb{R}\to \mathbb{R}^{2}$ is a trend
function,

- $s:\mathbb{R}\to \mathbb{R}$ is a scale function, and

- $(x_i^{(0)},y_i^{(0)})^T$ are elliptically distributed with zero location
$\mu_0 = 0$, scatter $\Sigma^{(0)}$ with $\mathrm{det}(\Sigma^{(0)}) = 1$, and
generating variate $\mathcal{R}^{(0)}$ with $\mathbb{E}(\mathcal{R}^{(0)}) = 1$.

We assume that both the trend and the scale depend on the day of the year. That
is, for each year the estimated trend and scale are the same. The trend and the
scale are estimated nonparametrically, and the estimates depend on the chosen
window size for the pooled rolling averages. Scatter plot of the standardized
coordinates $(x_i^{(0)}, y_i^{(0)})^T$, and the estimated trend and scale are
plotted. 

### Analysis of the data

The script `compute-k-and-acf.R` estimates the extreme value index of the
generating variate for various values of $k$. Furthermore, autocorrelation
function (acf) is computed. Original data has short-term dependence, and thus,
as a robustness check estimation and acf steps are repeated for thinned data
(every seventh day).

The script `compute-region.R` estimates extreme quantile regions for days $t$
with minimum and maximum estimates of $s(t)$ with chosen parameters.

### Running the empirical example

All the steps of the empirical example with chosen parameters are performed by
running the following.

```
bash main.sh
```

## Simulation study

Goal is to estimate $(1-1/n)$-quantile of the generating variate $\mathcal{R}$ of
a certain elliptical distribution. Simulations compare performance of four
different quantile estimators:

1. Extreme quantile estimator computed with

    i) true observations $\mathcal{R}_i$, and

    ii) approximations of the generating variate $\sqrt{(X_i -
    \hat\mu)^T\hat\Sigma_i(X_i - \hat\mu)}$.

2. Nonparametric estimator (largest order statistic) computed with true
   observations and approximations.

Boxplots of standardized estimation errors are outputted. Simulations with the
chosen parameters can be performed by running the following.
```
Rscript simulate.R
```
