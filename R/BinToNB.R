#' This function implements Step 5 of the algorithm.
#' It converts the multivariate binary data back to the original GPD scale
#'
#' @param prop.vec.bin vector of binary probabilities
#' @param NBprop NB proportion
#' @param Mlocation locations of median in the vector
#' @param bin.data generated multivariate binary data
#'
#' @return multivariate NB data and its correlation matrix
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
