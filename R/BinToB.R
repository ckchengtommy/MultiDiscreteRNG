#' Convert multivariate binary data back to the original binomial scale
#'
#' This function maps multivariate binary data to multivariate binomial outcomes,
#' preserving the original marginal distribution characteristics. Given the binary
#' representation of the data, the function  assigns binomial values
#' based on the original probability mass functions and the location of the median split.
#'
#' @importFrom stats cor
#' @param prop.vec.bin A vector of binary probabilities
#' @param BProp Binary proportion
#' @param Mlocation Indices of the medians in the vector
#' @param bin.data Generated multivariate binary data.
#' @return A list containing the multivariate binomial data and its correlation matrix
#' @examples
#' # Generate binary probabilities and probability mass functions for 3 variables
#' B.n.vec <- c(3, 4, 5)
#' B.prob.vec <- c(0.5, 0.5, 0.5)
#' p <- calc.bin.prob.B(B.n.vec, B.prob.vec)
#' pvec <- p$p
#' prop <- p$prop
#' Mlocation <- p$Mlocation
#'
#' # Select the first two variables for demonstration
#' pvec.pair      <- pvec[1:2]
#' Mlocation.pair <- Mlocation[1:2]
#' prop.pair      <- list(prop[[1]], prop[[2]])
#'
#' # Specify a target correlation matrix for two binary variables
#' del.next <- matrix(c(1.0, 0.3,
#'                      0.3, 1.0),
#'                    nrow = 2, byrow = TRUE)
#'
#' # Simulate N = 100 binary observations with the desired correlation
#' inter_bin <- generate.binaryVar(100, pvec.pair, del.next)
#'
#' # Convert back to binomial scale
#' Mydata <- BinToB(pvec.pair, prop.pair, Mlocation.pair, inter_bin)
#'
#' @export
#'
BinToB <- function(prop.vec.bin, BProp, Mlocation, bin.data){
  J = length(prop.vec.bin)
  K = numeric(0)
  prop = BProp

  B.data = matrix(NA, nrow(bin.data), ncol(bin.data))
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

    B.data[w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) < Mlocation[j]])), sum(w1),
                           prob = p[which(as.numeric(names(prop[[j]]))< Mlocation[j])], replace = TRUE)

    if(length(p[which(as.numeric(names(prop[[j]])) >=Mlocation[j])]) == 1){
      B.data[!w1, j] = as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) >= Mlocation[j]]))
    }else {
      B.data[!w1, j] = sample(as.numeric(names(prop[[j]][as.numeric(names(prop[[j]])) >= Mlocation[j]])), sum(!w1),
                              prob = p[which(as.numeric(names(prop[[j]])) >=Mlocation[j])], replace = TRUE)
    }
  }
  return(list(y = B.data, Corr = cor(B.data)))
}
