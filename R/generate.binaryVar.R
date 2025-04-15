#' Generate multivariate Binary data using the Emrich and Piedmonte (1991) approach
#' Approach
#'
#' @importFrom stats qnorm
#' @importFrom mvtnorm rmvnorm
#' @importFrom MultiOrd validation.CorrMat
#' @importFrom GenOrd corrcheck contord
#' @param nObs number of observations
#' @param prop.vec.bin probability of binary variables in a vector
#' @param corr.mat Correlation matrix
#' @return multivariate Binary Data
#' @export

generate.binaryVar <- function(nObs, prop.vec.bin, corr.mat)
{
  validation.CorrMat(prop.vec.bin, corr.mat)
  sigma_star = discrete_cont(as.list(prop.vec.bin), corr.mat)$SigmaC
  d = ncol(sigma_star)
  xx1 = rmvnorm(nObs, mean = rep(0, d), sigma = sigma_star)
  p = prop.vec.bin
  q = 1 - p
  BB= matrix(0, nObs, d)
  for (j in 1:d){
    for(i in 1:nObs){
      if (1 * xx1[i, j] > qnorm(1 - prop.vec.bin[j]))
        BB[i, j] = 1
      else BB[i, j] = 0
    }
  }
  return(BB)
}
