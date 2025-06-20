#' Convert Multivariate Binary Data Back to the Original Binomial Scale
#'
#' This function implements step 5 of the algorithm described in the paper.
#' It converts multivariate binary data back to the original binomial scale.
#'
#' @importFrom stats cor
#' @param prop.vec.bin A vector of binary probabilities
#' @param BProp Binary proportion
#' @param Mlocation Indices of the medians in the vector
#' @param bin.data Generated multivariate binary data.
#' @return A list containing the multivariate binomial data and its correlation matrix
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
