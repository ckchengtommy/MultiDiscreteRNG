#' Validate if the input GPD parameters are within feasible range
#'
#' This function returns the sum of two numbers.
#'
#' @param GPD.theta.vec Vector of theta values
#' @param GPD.lambda.vec Vector of lambda values
#' @return No return values; called it to check parameter inputs
#' @examples
#' validation.GPDparameters(GPD.theta.vec = c(3, 2), GPD.lambda.vec = c(0.4, 0.2))
#' @export

validation.GPDparameters <- function(GPD.theta.vec, GPD.lambda.vec)
{
  w1 = GPD.theta.vec <0
  w2 = GPD.lambda.vec >=1
  w3 = GPD.lambda.vec < 0 & GPD.lambda.vec < (-GPD.theta.vec)/4

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
