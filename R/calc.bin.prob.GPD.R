#' This function implements Step 1 of the algorithm.
#' It collapses the discrete outcome to binary ones for each variable.
#'
#' @param theta.vec vector of theta values
#' @param lambda.vec vector of lambda values
#'
#' @return vector of binary probability, dichotomous threshold
#'
#' @export
calc.bin.prob.GPD <- function(theta.vec, lambda.vec){

  validation.GPDparameters(theta.vec, lambda.vec)
  J <- length(theta.vec)
  p <- numeric(J)
  prop <- list(mode = vector)
  Mlocation <- numeric(J)

  for(i in 1:J){
    X.prop = GetGpoisPMF(1, theta.vec[i], lambda.vec[i])

    X_median = QuantileGpois(0.5, theta.vec[i], lambda.vec[i])

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
