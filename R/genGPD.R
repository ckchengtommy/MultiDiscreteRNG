#' Generate multivariate generalized Poisson data
#'
#' This function is the engine for simulating correlated generalized Poisson outcomes once the
#' intermediate binary parameters have been computed. It first generates correlated
#' multivariate binary data using a latent normal approach
#' implemented in \code{generate.binaryVar} and then
#' maps the binary outcomes back to the original generalized Poisson scales via a reverse-collapsing
#' step implemented in \code{BinToGPD}, using the generalized Poisson probability mass function and
#' dichotomization locations stored in \code{binObj}.
#'
#' @param no.rows integer; number of observations to generate (sample size \eqn{N}).
#' @param binObj list; intermediate object produced by the binary-correlation calibration
#'   step for generalized Poisson margins (e.g., \code{simBinaryCorr.GPD}). It must contain:
#'   \describe{
#'     \item{\code{pvec}}{numeric vector of binary success probabilities \eqn{p_j^b}.}
#'     \item{\code{intermat}}{intermediate/tetrachoric correlation matrix for the latent normal
#'       model used to generate correlated binary data.}
#'     \item{\code{GPDprop}}{list of generalized Poisson PMFs (probability masses over the support)
#'       used in the reverse-collapsing step.}
#'     \item{\code{Mlocation}}{numeric vector of dichotomization thresholds (median locations)
#'       used to split each marginal distribution into binary categories.}
#'   }
#'
#' @return A generated multivariate generalized Poisson dataset (returned as a list containing
#'   the simulated data matrix and its empirical correlation matrix).
#'
#' @examples
#' lambda.vec <- c(0.1, 0.2)
#' theta.vec <- c(7, 3)
#' M<- c(0.3, 0.3)
#' N <- diag(2)
#' N[lower.tri(N)] <- M
#' cmat<- N + t(N)
#' diag(cmat) <- 1
#'
#' # In real data simulation, no.rows should set to 100000 for accurate data generation
#' # in the intermediate step.
#' binObj = simBinaryCorr.GPD(GPD.theta.vec = theta.vec, GPD.lambda.vec = lambda.vec,
#'                            CorrMat = cmat, no.rows = 200, steps= 0.025)
#' data = genGPD(no.rows = 100, binObj = binObj)$y
#'
#' @export

genGPD <-function (no.rows, binObj)
{
  #browser()
  inter_bin = generate.binaryVar(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToGPD(binObj$pvec, binObj$GPDprop, binObj$Mlocation,
                    inter_bin)
  return(Mydata)
}
