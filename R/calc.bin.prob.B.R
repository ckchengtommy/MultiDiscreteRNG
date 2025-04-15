#' This function implements Step 1 of the algorithm.
#' It collapses the discrete outcome to binary ones for each variable.
#'
#' @param n.vec vector of number of trials
#' @param p.vec vector of probabilities
#'
#' @return vector of binary probability, dichotomization threshold
#'
#' @export


calc.bin.prob.B <- function(n.vec, p.vec){

  validation.Bparameters(n.vec, p.vec)
  J <- length(n.vec)
  p <- numeric(J)
  prop <- list(mode = vector)
  Mlocation <- numeric(J)

  for(i in 1:J){
    b.upper.limit = n.vec[i]
    X.prop = dbinom(c(0:b.upper.limit) , size = n.vec[i], prob = p.vec[i])
    names(X.prop) = as.character (c(0:b.upper.limit))
    X_median = qbinom(0.5, size = n.vec[i], prob = p.vec[i])
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
