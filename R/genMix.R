#' Generate multivariate mixed discrete data
#'
#' This function is the engine for simulating correlated mixed discrete outcomes once the
#' intermediate binary parameters have been computed. It first generates correlated
#' multivariate binary data using a latent normal approach
#' implemented in \code{generate.binaryVar} and then
#' maps the binary outcomes back to the original mixed discrete scales via a reverse-collapsing
#' step implemented in \code{BinToMix}, using the component probability mass functions and
#' dichotomization locations stored in \code{binObj}.
#'
#' @param no.rows integer; number of observations to generate (sample size \eqn{N}).
#' @param binObj list; intermediate object produced by the binary-correlation calibration
#'   step for mixed discrete margins (e.g., \code{simBinaryCorr.Mix}). It must contain:
#'   \describe{
#'     \item{\code{pvec}}{numeric vector of binary success probabilities.}
#'     \item{\code{intermat}}{intermediate/tetrachoric correlation matrix for the latent normal
#'       model used to generate correlated binary data.}
#'     \item{\code{Mixprop}}{list of probability mass functions (probability masses over the support)
#'       for each discrete margin, used in the reverse-collapsing step.}
#'     \item{\code{Mlocation}}{numeric vector of dichotomization thresholds (median locations)
#'       used to split each marginal distribution into binary categories.}
#'   }
#'
#' @return A generated multivariate mixed discrete dataset (returned as a list containing
#'   the simulated data matrix and its empirical correlation matrix).
#'
#' @examples
#' #Define parameters
#' GPD.theta <- 2
#' GPD.lambda <- 0.3
#' NB.r <- 10
#' NB.prob <- 0.2
#' B.n <- 5
#' B.prob <- 0.5
#'
#'
#' #Define Correlation Matrix
#' M<- c(0.3, 0.3, 0.3)
#' N <- diag(3)
#' N[lower.tri(N)] <- M
#' cmat<- N + t(N)
#' diag(cmat) <- 1
#'
#'
#' # In real data simulation, no.rows should set to 100000 for accurate data generation
#' # in the intermediate step.
#' binObj <- simBinaryCorr.Mix(GPD.theta.vec = GPD.theta, GPD.lambda.vec = GPD.lambda,
#'   NB.r.vec = NB.r, NB.prob.vec = NB.prob,
#'   B.n.vec = B.n, B.prob.vec = B.prob,
#'   CorrMat = cmat, no.rows = 250)
#'
#' MixData = genMix(no.rows = 100, binObj = binObj)$y
#'
#' @export


genMix = function (no.rows, binObj)
{
  #browser()
  inter_bin = generate.binaryVar(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToMix(binObj$pvec, binObj$Mixprop, binObj$Mlocation,
                    inter_bin)
  return(Mydata)
}
