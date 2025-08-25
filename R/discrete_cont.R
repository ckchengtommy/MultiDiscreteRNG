#' Compute the tetrachoric correlation matrix for a multivariate standard normal distribution
#'
#' This function calculates the intermediate correlation matrix of a multivariate
#' standard normal distribution in Step 2 of the algorithm. If the resulting matrix
#' is not positive definite, the nearest positive definite matrix is returned and
#' a warning is issued.
#'
#' @importFrom Matrix nearPD
#' @importFrom GenOrd corrcheck
#' @param marginal a list of \eqn{k} elements, where \eqn{k} is the number of variables.
#' The \eqn{i}-th element of \code{marginal} is the vector of the cumulative probabilities defining the marginal distribution of the \eqn{i}-th component of the  multivariate variable. If the \eqn{i}-th component can take \eqn{k_i} values, the \eqn{i}-th element of \code{marginal} will contain \eqn{k_i-1} probabilities (the \eqn{k_i}-th is obviously 1 and shall not be included).
#' @param Sigma the target correlation matrix of the discrete variables
#' @param Spearman A logical flag indicating whether Spearman correlation should be used
#' @param maxit maximum iterations of the algorithm to correct PD matrix
#' @param epsilon tolerance of the algorithm convergence
#' @param support a list of \eqn{k} elements, where \eqn{k} is the number of variables. The \eqn{i}-th element of \code{support} is the vector containing the ordered values of the support of the \eqn{i}-th variable. By default, the support of the \eqn{i}-th variable is \eqn{1,2,...,k_i}
#' @references Ferrari and Barbiero 2012 (<https://doi.org/10.1080/00273171.2012.692630>)
#' @return No return values; called it to check parameter inputs
#' @examples
#' prop.vec.bin = c(0.5037236, 0.5034147)
#' cor.mat = matrix(c(1, 0.3, 0.3, 1), nrow = 2, byrow = TRUE)
#' InterMVN_Sigma = discrete_cont(marginal = prop.vec.bin, Sigma = cor.mat)$SigmaC
#'
#' @export


discrete_cont <- function (marginal, Sigma, support = list(), Spearman = FALSE,
                            epsilon = 1e-06, maxit = 100)
{
  if (!all(unlist(lapply(marginal, function(x) (sort(x) ==
                                                x & min(x) > 0 & max(x) < 1)))))
    stop("Error in assigning marginal distributions!")
  if (!isSymmetric(Sigma) | min(eigen(Sigma)$values) < 0 |
      !all(diag(Sigma) == 1))
    stop("Correlation matrix not valid!")
  len <- length(support)
  k <- length(marginal)
  niter <- matrix(0, k, k)
  kj <- numeric(k)
  for (i in 1:k) {
    kj[i] <- length(marginal[[i]]) + 1
    if (len == 0) {
      support[[i]] <- 1:kj[i]
    }
  }
  for (g in 2:k) {
    if (det(Sigma[1:g, 1:g]) <= 0) {
      stop("Main minor number ", g, " is not positive!",
           "\n")
    }
  }
  mcmin <- corrcheck(marginal, support, Spearman)[[1]]
  mcmax <- corrcheck(marginal, support, Spearman)[[2]]
  if (sum(mcmin <= Sigma & Sigma <= mcmax) != k^2) {
    stop("Some correlation coefficients are not feasible!",
         "\n", "Please use function corrcheck to get lower and upper bounds!",
         "\n")
  }
  Sigma0 <- Sigma
  Sigmaord <- Sigma
  Sigmaord <- contord(marginal, Sigma, support, Spearman)
  Sigmaold <- Sigma
  Sigmaordold <- Sigmaord
  for (q in 1:(k - 1)) {
    for (r in (q + 1):k) {
      if (Sigma0[q, r] == 0) {
        Sigma[q, r] <- 0
      }
      else {
        it <- 0
        while (max(abs(Sigmaord[q, r] - Sigma0[q, r]) >
                   epsilon) & it < maxit) {
          if (Sigma0[q, r] * (Sigma0[q, r]/Sigmaordold[q,
                                                       r]) >= 1) {
            Sigma[q, r] <- Sigmaold[q, r] * (1 + 0.1 *
                                               (1 - Sigmaold[q, r]) * sign(Sigma0[q, r] -
                                                                             Sigmaord[q, r]))
          }
          else {
            Sigma[q, r] <- Sigmaold[q, r] * (Sigma0[q,
                                                    r]/Sigmaord[q, r])
          }
          Sigma[r, q] <- Sigma[q, r]
          Sigmaord[r, q] <- contord(list(marginal[[q]],
                                         marginal[[r]]), matrix(c(1, Sigma[q, r],
                                                                  Sigma[q, r], 1), 2, 2), list(support[[q]],
                                                                                               support[[r]]), Spearman)[2]
          Sigmaord[q, r] <- Sigmaord[r, q]
          Sigmaold[q, r] <- Sigma[q, r]
          Sigmaold[r, q] <- Sigmaold[q, r]
          it <- it + 1
        }
        niter[q, r] <- it
        niter[r, q] <- it
      }
    }
  }
  if (eigen(Sigma)$values[k] <= 0) {
    warning("The Sigma matrix is not coherent with the given margins or\nit is impossible to find a feasible correlation matrix for MVN\nensuring Sigma for the given margins")
    Sigma = as.matrix(nearPD(Sigma, corr = TRUE, conv.norm.type = "I", keepDiag = TRUE)$mat)
    Sigma = (Sigma + t(Sigma))/2
  }
  emax <- max(abs(Sigmaord - Sigma0))
  list(SigmaC = Sigma, SigmaO = Sigmaord, Sigma = Sigma0, niter = niter,
       maxerr = emax)
}
