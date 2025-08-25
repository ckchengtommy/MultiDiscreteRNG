#' Validate if the input NB parameters are within feasible range
#'
#'
#' @param NB.r.vec Vector of number of trials parameters
#' @param NB.prob.vec Vector of probabilities
#' @return No return values; called it to check parameter inputs
#' @examples
#' validation.NBparameters(NB.r.vec = c(10, 15), NB.prob.vec = c(0.7, 0.5))
#' @export


validation.NBparameters <- function(NB.r.vec, NB.prob.vec)
{
  w1 = NB.r.vec <0
  w2 = NB.prob.vec <0
  w3 = NB.prob.vec >1

  if(sum(w1)>0){
    stop("number of trials must be greater than 0")
  }
  if(sum(w2)>0){
    stop("The probabilities cannot be less than 0")
  }
  if(sum(w3)>0){
    stop("The probabilities cannot be greater than 1")
  }
}
