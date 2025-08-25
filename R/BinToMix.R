#' Convert multivariate binary data to mixed distribution outcomes
#'
#' This function maps multivariate binary data to outcomes from a mixture of
#' generalized Poisson, negative binomial, and binomial distributions while
#' preserving the original marginal distribution characteristics. It assigns
#' appropriate values from each distribution based on binary thresholds and
#' probability mass functions.
#'
#' @param prop.vec.bin A vector of binary probabilities for each variable
#' @param Mixprop A list of probability mass functions for each variable's distribution
#' @param Mlocation A vector of threshold values (typically medians) for each variable
#' @param bin.data A matrix of multivariate binary data (0s and 1s)
#'
#' @return A list containing:
#' \item{y}{Matrix of generated mixed distribution data}
#' \item{Corr}{Correlation matrix of the generated data}
#'
#' @examples
#' \dontrun{
#' # First simulate intermediate binary correlations
#' result <- simBinaryCorr.Mix(
#' GPD.theta.vec = c(1, 2),
#' GPD.lambda.vec = c(0.5, 0.3),
#' NB.r.vec = 10,
#' NB.prob.vec = 0.2,
#' B.n.vec = 5,
#' B.prob.vec = 0.5,
#' CorrMat = matrix(c(1, 0.3, 0.2, 0.1,
#' 0.3, 1, 0.4, 0.2,
#' 0.2, 0.4, 1, 0.3,
#' 0.1, 0.2, 0.3, 1), 4, 4),
#' no.rows = 1000
#' )
#'
#' # Generate correlated binary data using the intermediate matrix
#' bin_data <- generate.binary(1000, result$pvec, result$intermat)
#'
#' # Convert binary data to mixed distribution outcomes
#' mixed_data <- BinToMix(result$pvec, result$Mixprop, result$Mlocation, bin_data)
#' }
#'
#' @export


BinToMix = function(prop.vec.bin, Mixprop, Mlocation, bin.data){
  #browser()
  J = length(prop.vec.bin)
  K = numeric(0)
  prop = Mixprop

  Mix.data = matrix(NA, nrow(bin.data), ncol(bin.data))
  for(a in 1:length(prop)){
    K[a] = length(prop[[a]])
  }
  for (j in 1:J) {
    p = numeric(0)
    for (k in 1:K[j]) {
      if (k < which(names(prop[[j]]) == Mlocation[j])) {
        p[k] = prop[[j]][k]/(1-prop.vec.bin[j])
      }
      else {
        p[k] = prop[[j]][k]/ (1-prop.vec.bin[j])
      }
    }

    w1 = bin.data[, j] == 0

    Mix.data[w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) < Mlocation[j]])), sum(w1),
                             prob = p[which(as.numeric(names(prop[[j]]))< Mlocation[j])], replace = TRUE)

    Mix.data[!w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) >= Mlocation[j]])), sum(!w1),
                              prob = p[which(as.numeric(names(prop[[j]])) >=Mlocation[j])], replace = TRUE)
  }
  return(list(y = Mix.data, Corr = cor(Mix.data)))
}
