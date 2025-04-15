#' Validate if the input Binomial parameters are within feasible range
#'
#' This function returns the sum of two numbers.
#'
#' @param n.vec Vector of number of trials
#' @param p.vec Vector of probability
#' @return No return values; called it to check parameter inputs
#' @examples
#' validation.Bparameters(n.vec = c(10, 15), p.vec = c(0.4, 0.2))
#' @export


validation.Bparameters <- function(n.vec, p.vec)
{
  n.vec <- floor(n.vec)
  w1 = n.vec <0
  w2 = p.vec <0
  w3 = p.vec >1

  if(sum(w1)>0){
    stop("Number of trials must be greater than 0")
  }
  if(sum(w2)>0){
    stop("The probabilities cannot be less than 0")
  }
  if(sum(w3)>0){
    stop("The probabilities cannot be greater than 1")
  }
}
