#' This function implements Step 5 of the algorithm
#' It converts the multivariate binary data back to the original GPD scale
#'
#' @param prop.vec.bin vector of binary probabilities
#' @param GPDprop GPD proportion
#' @param Mlocation locations of median in the vector
#' @param bin.data generated multivariate binary data
#'
#' @return multivariate GPD data and its correlation matrix
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
