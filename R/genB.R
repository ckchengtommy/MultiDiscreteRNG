#' Generate multivariate binomial data
#'
#' This function is the engine for simulating correlated binomial outcomes once the
#' intermediate binary parameters have been computed. It first generates correlated
#' multivariate binary data using a latent normal approach implemented in \code{generate.binaryVar} and then
#' maps the binary outcomes back to the original binomial scales via a reverse-collapsing
#' step implemented in \code{BinToB}, using the binomial probability mass function and dichotomization
#' locations stored in \code{binObj}
#'
#' @param no.rows integer; number of observations to generate (sample size \eqn{N}).
#' @param binObj list; intermediate object produced by the binary-correlation calibration
#'   step for binomial margins (e.g., \code{simBinaryCorr.B}). It must contain:
#'   \describe{
#'     \item{\code{pvec}}{numeric vector of binary success probabilities \eqn{p_j^b}.}
#'     \item{\code{intermat}}{intermediate/tetrachoric correlation matrix for the latent normal
#'       model used to generate correlated binary data.}
#'     \item{\code{BProp}}{list of binomial PMFs (probability masses over the support) used in
#'       the reverse-collapsing step.}
#'     \item{\code{Mlocation}}{numeric vector of dichotomization thresholds (median locations)
#'       used to split each marginal distribution into binary categories.}
#'   }
#'
#' @return A generated multivariate binomial dataset (returned as a list containing
#'   the simulated data matrix and its empirical correlation matrix.
#'
#' @examples
#' n.vec <- c(3, 4, 5)
#' p.vec <- c(0.5, 0.5, 0.5)
#'
#' M<- c(0.3, 0.4, 0.3)
#' N <- diag(3)
#' N[lower.tri(N)] <- M
#' cmat<- N + t(N)
#' diag(cmat) <- 1
#' #In real data simulation, no.rows should set to 100000 for accurate data generation
#' #in the intermediate step
#' binObj = simBinaryCorr.B(B.n.vec = n.vec, B.prob.vec = p.vec,
#' CorrMat = cmat, no.rows = 20000, steps= 0.025)
#'
#' data = genB(no.rows = 100, binObj = binObj)$y
#' @export

genB <-function (no.rows, binObj)
{
  inter_bin = generate.binaryVar(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToB(binObj$pvec, binObj$BProp, binObj$Mlocation,
                  inter_bin)
  return(Mydata)
}
