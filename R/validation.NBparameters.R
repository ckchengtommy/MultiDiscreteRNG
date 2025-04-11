#' Validate if the input NB parameters are within feasible range
#'
#'
#' @param r.vec Vector of number of trials parameters
#' @param prob.vec Vector of probabilties
#' @return No return values; called it to check parameter inputs
#' @examples
#' validation.NBparameters(r.vec = c(10, 15), prob.vec = c(0.7, 0.5))
#' @export


validation.NBparameters <- function(r.vec, prob.vec)
{
  w1 = r.vec <0
  w2 = prob.vec <0
  w3 = prob.vec >1

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
