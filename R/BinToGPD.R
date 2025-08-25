#' Convert multivariate binary data back to the original generalized Poisson scale
#'
#' This function maps multivariate binary data to multivariate generalized Poisson outcomes,
#' preserving the original marginal distribution characteristics. Given a binary representation,
#' it assigns generalized Poisson values based on the original probability mass
#' functions and the location of the median split for each variable.
#'
#' @param prop.vec.bin A vector of binary probabilities
#' @param GPDprop Generalized Poisson distribution probability mass functions tables
#' @param Mlocation Indices of the medians in the vector
#' @param bin.data Generated multivariate binary data
#'
#' @return A list containing the multivariate generalized Poisson data and its correlation matrix
#' @examples
#' # Prepare three GPD parameter vectors
#' GPD.lambda.vec <- c(0.1, 0.2, 0.3)
#' GPD.theta.vec  <- c(7, 0.7, 40)
#'
#' # Compute binary probabilities, PMFs, and thresholds
#' p     <- calc.bin.prob.GPD(GPD.theta.vec, GPD.lambda.vec)
#' pvec  <- p$p
#' prop  <- p$prop
#' Mloc  <- p$Mlocation
#'
#' # Use only the first two variables for demonstration
#' pvec.pair      <- pvec[1:2]
#' Mlocation.pair <- Mloc[1:2]
#' prop.pair      <- list(prop[[1]], prop[[2]])
#'
#' # Define a 2×2 target correlation matrix
#' del.next <- matrix(c(1.0, 0.3,
#'                      0.3, 1.0),
#'                    nrow = 2, byrow = TRUE)
#'
#' # Simulate 100 correlated binary observations
#' inter_bin <- generate.binaryVar(100, pvec.pair, del.next)
#'
#' # Reconstruct the GPD‐scaled data
#' Mydata <- BinToGPD(pvec.pair, prop.pair, Mlocation.pair, inter_bin)
#'
#' @export
BinToGPD <- function(prop.vec.bin, GPDprop, Mlocation, bin.data){
  J = length(prop.vec.bin)
  K = numeric(0)
  prop = GPDprop

  GPD.data = matrix(NA, nrow(bin.data), ncol(bin.data))
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
        p[k] = prop[[j]][k]/ (1-prop.vec.bin[j])
      }
    }

    w1 = bin.data[, j] == 0

    GPD.data[w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) < Mlocation[j]])), sum(w1),
                             prob = p[which(as.numeric(names(prop[[j]]))< Mlocation[j])], replace = TRUE)

    GPD.data[!w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) >= Mlocation[j]])), sum(!w1),
                              prob = p[which(as.numeric(names(prop[[j]])) >=Mlocation[j])], replace = TRUE)
  }
  return(list(y = GPD.data, Corr = cor(GPD.data)))
}
