#' Validate binomial parameters are within feasible ranges
#'
#' This helper function checks that binomial inputs are valid before downstream
#' probability calculations and data generation. If any check fails, the function stops
#' with an informative error message.
#'
#' @param B.n.vec Vector of number of trials
#' @param B.prob.vec Vector of probability
#' @return No return values; called it to check parameter inputs
#' @examples
#' validation.Bparameters(B.n.vec = c(10, 15), B.prob.vec = c(0.4, 0.2))
#' @export


validation.Bparameters <- function(B.n.vec, B.prob.vec)
{
  B.n.vec <- floor(B.n.vec)
  w1 = B.n.vec <0
  w2 = B.prob.vec <0
  w3 = B.prob.vec >1

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
