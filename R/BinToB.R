#' This function implements Step 5 of the algorithm
#' It converts the multivariate binary data back to the original Binomial scale
#'
#' @param prop.vec.bin vector of binary probabilities
#' @param Bprop B proportion
#' @param Mlocation locations of median in the vector
#' @param bin.data generated multivariate binary data
#'
#' @return multivariate Binomial data and its correlation matrix
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
