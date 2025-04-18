## Overview

The R package MultiDiscreteRNG implements the binary collapsing and re-collapsing approach by Demirtas 2006 to generate Multivariate correlated discrete data with generalized Poisson, negative binomial and binomial marginal distributions.

## Installation

``` r
devtools::install_github("ckchengtommy/MultiDiscreteRNG", force = TRUE)
library(MultiDiscreteRNG)
```

You can also install the official release version from CRAN

``` r
install.packages("MultiDiscreteRNG")
```

## Example

1.  Generate a trivariate generalized Poisson distribution with dispersion parameters $\lambda_1 = 0.1, \lambda_2 = 0.2, \lambda_3 = 0.3$, rate parameters $\theta_1 = 7, \theta_2 = 0.7, \theta_3 = 40$, with a correlation matrix, and pairwise correlation $\rho_{12} = 0.3, \rho_{13} = 0.3, \rho_{23} = 0.3$

First, specify the marginal distribution parameters with

``` r
lambda.vec <- c(0.1, 0.2, 0.3) 
theta.vec <- c(7, 0.7, 40)
```

Create the correlation matrix called cmat

``` r
M<- c(0.3, 0.3, 0.3)
N <- diag(3)
N[lower.tri(N)] <- M
cmat<- N + t(N)
diag(cmat) <- 1
```

Then we generate the intermediate correlation matrix with

``` r
binObj = simBinaryCorr.GPD(theta.vec, lambda.vec, cmat, 1000000, steps= 0.025)
```

Lastly, generate the trivaiarte generalized Poisson distribution with n = 100 samples

``` r
GPD.data = genGPD(no.rows = 100, binObj)
```
