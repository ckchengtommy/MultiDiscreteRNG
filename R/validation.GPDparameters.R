#' Validate if the input GPD parameters are within feasible range
#'
#' This function returns the sum of two numbers.
#'
#' @param theta.vec Vector of theta values
#' @param lambda.vec Vector of lambda values
#' @return No return values; called it to check parameter inputs
#' @examples
#' validation.GPDparameters(theta.vec = c(3, 2), lambda.vec = c(0.4, 0.2))
#' @export

validation.GPDparameters <- function(theta.vec, lambda.vec)
{
  w1 = theta.vec <0
  w2 = lambda.vec >=1
  w3 = lambda.vec < 0 & lambda.vec < (-theta.vec)/4

  if(sum(w1)>0){
    stop("theta has to be greater than 0!")
  }
  if(sum(w2)>0){
    stop("lambda has to be less than 1!")
  }
  if(sum(w3)>0){
    stop(paste("For lambda < 0, lambda must be greater than or equal to -theta/4"))
  }
}
