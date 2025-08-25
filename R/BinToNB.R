#' Convert multivariate binary data back to the original negative binomial scale
#'
#' This function maps multivariate binary data to multivariate negative binomial outcomes,
#' preserving the original marginal distribution characteristics. Given a binary representation,
#' the function assigns negative binomial values based on the original probability
#' mass functions and the location of the median split for each variable.
#'
#' @param prop.vec.bin A numeric vector of binary probabilities
#' @param NBprop A numeric value or vector of negative binomial proportions
#' @param Mlocation Integer indices of the medians in the vector
#' @param bin.data A data frame or matrix of generated multivariate binary data
#'
#' @return A list containing the multivariate negative binomial data and its correlation matrix
#' @examples
#' NB.r.vec <- c(10, 3, 16)
#' NB.prob.vec <- c(0.65, 0.4, 0.88)
#'
#' # Compute binary probabilities, PMFs, and thresholds
#' p      <- calc.bin.prob.NB(NB.r.vec, NB.prob.vec)
#' pvec   <- p$p
#' prop   <- p$prop
#' Mloc   <- p$Mlocation
#'
#' # Use only the first two variables for demonstration
#' pvec.pair      <- pvec[1:2]
#' Mlocation.pair <- Mloc[1:2]
#' prop.pair      <- list(prop[[1]], prop[[2]])
#'
#' # Define a 2×2 target correlation matrix
#' del.next <- matrix(c(1.0, -0.3,
#'                      -0.3, 1.0),
#'                   nrow = 2, byrow = TRUE)
#'
#' # Simulate 100 correlated binary observations
#' inter_bin <- generate.binaryVar(100, pvec.pair, del.next)
#'
#' # Reconstruct the negative binomial scaled data
#' Mydata <- BinToNB(pvec.pair, prop.pair, Mlocation.pair, inter_bin)
#'
#' @export
BinToNB <- function(prop.vec.bin, NBprop, Mlocation, bin.data){
  J = length(prop.vec.bin)
  K = numeric(0)
  prop = NBprop

  NB.data = matrix(NA, nrow(bin.data), ncol(bin.data))
  for(a in 1:length(prop)){
    K[a] = length(prop[[a]])
  }
  for (j in 1:J) {
    p = numeric(0)
    for (k in 1:K[j]) {
      if (k < which(names(prop[[j]])==Mlocation[j])) {
        p[k] = prop[[j]][k]/(1-prop.vec.bin[j])
      }
      else {
        p[k] = prop[[j]][k]/ (prop.vec.bin[j])
      }
    }

    w1 = bin.data[, j] == 0

    NB.data[w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) < Mlocation[j]])), sum(w1),
                            prob = p[which(as.numeric(names(prop[[j]]))< Mlocation[j])], replace = TRUE)

    NB.data[!w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) >= Mlocation[j]])), sum(!w1),
                             prob = p[which(as.numeric(names(prop[[j]])) >=Mlocation[j])], replace = TRUE)
  }
  return(list(y = NB.data, Corr = cor(NB.data)))
}
