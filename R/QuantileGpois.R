#' This function computes the quantile of Generalized Poisson
#'
#'
#' @param p vector of probabilities
#' @param theta vector of theta
#' @param lambda vector of lambda
#'
#' @return the quantile of GPD
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
    cat(paste0("x = 0, P(X = x) = ", round(p, 7), ", P(X <= x) = ",
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
        cat(paste0("x = ", i, ", P(X = x) = ", round(p,
                                                     7), ", P(X <= x) = ", round(s[i + 1], 7), "\n"))
      i = i + 1
    }
  }
  if (lambda < 0) {
    Fm = s[i]
    s[1:i] = s[1:i] * Fm^(-1)
    if (details) {
      cat("When lambda is negative, we need to account for truncation error. ")
      cat("The adjusted CDF are:", s[1:i], "\n")
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
