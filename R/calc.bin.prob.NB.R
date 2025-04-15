#' This function implements Step 1 of the algorithm.
#' It collapses the discrete outcome to binary ones for each variable.
#'
#' @importFrom stats rnbinom dnbinom qnbinom
#' @param r.vec vector of number of trials
#' @param prob.vec vector of probabilities
#'
#' @return vector of binary probability, dichotomous threshold
#'
#' @export

calc.bin.prob.NB <- function(r.vec, prob.vec){
  validation.NBparameters(r.vec, prob.vec)
  J <- length(r.vec)
  p <- numeric(J)
  prop <- list(mode = vector)
  Mlocation <- numeric(J)

  for(i in 1:J){
    j = rnbinom(100000, size = r.vec[i], prob = prob.vec[i])
    j.table = table(j)/100000
    nb.upper.limit = as.numeric(names(j.table)[length(j.table)])

    while(  dnbinom(nb.upper.limit, size = r.vec[i], prob = prob.vec[i]) >= 1e-10){
      nb.upper.limit = nb.upper.limit +1
    }

    X.prop = dnbinom(c(0:nb.upper.limit) , size = r.vec[i], prob = prob.vec[i])
    names(X.prop) = as.character (c(0:nb.upper.limit))
    X_median = qnbinom(0.5, size = r.vec[i], prob = prob.vec[i])
    p0 = sum(X.prop[as.numeric(names(X.prop)) < X_median])
    p1 = sum(X.prop[as.numeric(names(X.prop)) > X_median])
    p_median = X.prop[as.numeric(names(X.prop)) == X_median]

    if( abs(p0+p_median -0.5) > abs(p1+p_median-0.5)){
      p1 <- p1+p_median
      Mlocation[i] <- X_median
    }else{
      p0 <- p0+p_median
      Mlocation[i] <- X_median+1
    }

    p[i] <- p1
    prop[[i]] <- X.prop
  }
  return(list(p = p, prop = prop, Mlocation = Mlocation))
}
