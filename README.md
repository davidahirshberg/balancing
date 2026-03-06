# balancing

> **Warning:** This package is in early development. The code was largely written by an AI assistant and has not been fully reviewed or tested by the author. Use at your own risk.

Kernel balancing weights via Bregman method-of-moments with pluggable dispersions, kernels, and estimands. Supports survival and single-outcome settings with adaptive regularization (Lepski) and cross-fitting.

## Installation

```r
# install.packages("remotes")
remotes::install_github("davidahirshberg/balancing")
```

## Usage

```r
library(balancing)

# Single-outcome ATE with entropy balancing
est <- balancing_estimate(
  Y = Y, A = A, X = X,
  estimand = treatment_specific_mean(),
  kernel = matern_kernel(nu = 3/2),
  dispersion = entropy_dispersion(),
  tuning = lepski_tuning()
)
```

## License

MIT
