#' Compute the quantile function of the generalized Poisson distribution
#'
#' This function evaluates the generalized Poisson quantile \eqn{Q(p)} by
#' incrementally constructing the PMF and CDF from \eqn{x=0} upward until the
#' CDF exceeds the largest requested probability in \code{p}. The returned
#' quantile(s) are the smallest integer \eqn{x} such that \eqn{P(X \le x) \ge p}.
#'
#' @param p vector of probabilities
#' @param theta vector of theta
#' @param lambda vector of lambda
#' @param details A logical flag to return the computational details
#' @return An integer vector of generalized Poisson quantiles corresponding to \code{p}.
#'
#' @examples
#' QuantileGpois(p = 0.95, theta = 2, lambda = 0.1)
#'
#' @export
#'
QuantileGpois <- function (p, theta, lambda, details = FALSE)
{
  if (theta <= 0)
    stop("theta has to be greater than 0!")
  if (lambda > 1)
    stop("lambda has to be less than 1!")
  if (lambda < 0 & lambda < (-theta)/4)
    stop(paste("For lambda < 0, lambda must be greater than or equal -theta/4, which is ",
               (-theta)/4, "!", sep = ""))
  if (max(p) > 1 | min(p) < 0)
    stop("p should be between 0 and 1!")
  m = numeric(1)
  if (lambda < 0) {
    mod = floor(-theta/lambda)
    if (-theta - mod * lambda == 0)
      m = mod
    else m = mod + 1
  }
  p.in = p
  upper = max(p.in)
  s = numeric(10000)
  q = numeric(length(p.in))
  w = exp(-lambda)
  p = exp(-theta)
  s[1] = p
  if (details)
    message(paste0("x = 0, P(X = x) = ", round(p, 7), ", P(X <= x) = ",
               round(s[1], 7), "\n"))
  i = 1
  while (s[i] < upper) {
    if (lambda < 0 & i > m)
      break
    else {
      p = theta * (theta + lambda * i)^(i - 1) * exp(-theta -
                                                       lambda * i)/factorial(i)
      if (i == 10000) {
        temp = numeric(10000)
        s = c(s, temp)
      }
      s[i + 1] = s[i] + p
      if (p == 0 | p == Inf)
        break
      if (details)
        message(paste0("x = ", i, ", P(X = x) = ", round(p,
                                                     7), ", P(X <= x) = ", round(s[i + 1], 7), "\n"))
      i = i + 1
    }
  }
  if (lambda < 0) {
    Fm = s[i]
    s[1:i] = s[1:i] * Fm^(-1)
    if (details) {
      message("When lambda is negative, we need to account for truncation error. ")
      message("The adjusted CDF are:", s[1:i], "\n")
    }
  }
  for (j in 1:length(p.in)) {
    i = 1
    while (p.in[j] > s[i]) {
      if (s[i] == 0 | s[i] == Inf)
        break
      i = i + 1
    }
    q[j] = i - 1
  }
  if (s[i] == 0 | s[i] == Inf)
    return(q - 1)
  else return(q)
}
