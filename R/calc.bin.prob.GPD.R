#' Collapse discrete generalized Poisson outcomes to binary variables
#'
#' This function implements Step 1 of the algorithm. It collapses each discrete
#' outcome from the generalized Poisson distribution (GPD) into a binary
#' probability and determines a dichotomous threshold for each variable.
#'
#' @param GPD.theta.vec Numeric vector of GPD theta parameters.
#' @param GPD.lambda.vec Numeric vector of GPD lambda parameters.
#' @return A list containing the binary probability vector, the list of GPD probability mass functions
#'   for each variable, and the corresponding threshold indices.
#' @examples
#' # Prepare three GPD parameter vectors
#' GPD.lambda.vec <- c(0.1, 0.2, 0.3)
#' GPD.theta.vec  <- c(7, 0.7, 40)
#'
#' # Compute binary probabilities, PMFs, and thresholds
#' p <- calc.bin.prob.GPD(GPD.theta.vec, GPD.lambda.vec)
#'
#'
#' @export
calc.bin.prob.GPD <- function(GPD.theta.vec, GPD.lambda.vec){

  validation.GPDparameters(GPD.theta.vec, GPD.lambda.vec)
  J <- length(GPD.theta.vec)
  p <- numeric(J)
  prop <- list(mode = vector)
  Mlocation <- numeric(J)

  for(i in 1:J){
    X.prop = GetGpoisPMF(1, GPD.theta.vec[i], GPD.lambda.vec[i])

    X_median = QuantileGpois(0.5, GPD.theta.vec[i], GPD.lambda.vec[i])

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
