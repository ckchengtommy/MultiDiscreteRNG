## Overview

The R package MultiDiscreteRNG implements the binary collapsing and re-collapsing approach by Demirtas 2006 to generate multivariate correlated discrete data with generalized Poisson, negative binomial and binomial marginal distributions.

## Installation

``` r
devtools::install_github("ckchengtommy/MultiDiscreteRNG", force = TRUE)
library(MultiDiscreteRNG)
```

You can also install the official release version from CRAN

``` r
install.packages("MultiDiscreteRNG")
```

## Example 1: Trivariate Generalized Poisson Distribution

To generate a trivariate generalized Poisson distribution with dispersion parameters $\lambda_1 = 0.1, \lambda_2 = 0.2, \lambda_3 = 0.3$, rate parameters $\theta_1 = 7, \theta_2 = 0.7, \theta_3 = 40$, with a correlation matrix $\rho_{12} = 0.3, \rho_{13} = 0.3, \rho_{23} = 0.3$

First, specify the marginal distribution parameters with

``` r
lambda.vec <- c(0.1, 0.2, 0.3) 
theta.vec <- c(7, 0.7, 40)
```

Create the correlation matrix called cmat:

``` r
M<- c(0.3, 0.3, 0.3)
N <- diag(3)
N[lower.tri(N)] <- M
cmat<- N + t(N)
diag(cmat) <- 1
```

Then we generate the intermediate correlation matrix with:

``` r
binObj = simBinaryCorr.GPD(theta.vec, lambda.vec, cmat, 1000000, steps= 0.025)
```

Lastly, generate the trivaiarte generalized Poisson distribution with n = 100 samples.

``` r
GPD.data = genGPD(no.rows = 100, binObj)$y
```

## Example 2: 4-variate Negative Binomial Distribution

Note: Let $y$ be the marginal negative binomial random variable. We define the parameterization of the negative binomial distribution as follows: y is the number of failures before $r^{th}$ success, n is the number of trials per observation and p is the probability of success.

To generate a 4-variate negative binomial distribution with number of trials parameter $r_1 = 3, r_2 = 6, r_3 = 9, r_4 = 10$, probability of success parameters $p_1 = 0.51, p_2 = 0.32, p_3 = 0.64, p_4 = 0.43$ and a unstructured correlation matrix with some arbitrary values.

First, specify the marginal distribution parameters with:

``` r
r.vec <- c(3, 6, 9, 10)
p.vec <- c(0.51, 0.32, 0.64, 0.43)
```

Create a correlation matrix with some arbitrary values:

``` r
M<- c(-0.0544 ,0.0381, -0.2319 , 0.0145 , 0.0758 ,0.1048)
N <- diag(4)
N[lower.tri(N)] <- M
cmat<- N + t(N)
diag(cmat) <- 1
```

Then we generate the intermediate correlation matrix with:

``` r
binObj = simBinaryCorr.NB(r.vec, p.vec, cmat, 1000000, steps= 0.025)
```

Lastly, generate the 4-variate negative binomial distribution with n = 100 samples:

``` r
NB.data = genNB(no.rows = 100, binObj)$y
```

A multivariate binomial distribution can be generated using a similar procedure.
