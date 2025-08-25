#' Collapse discrete negative binomial outcomes to binary variables
#'
#' This function implements Step 1 of the algorithm. It collapses each discrete
#' outcome from the negative binomial distribution into a binary
#' probability and determines a dichotomous threshold for each variable.
#'
#' @importFrom stats rnbinom dnbinom qnbinom
#' @param NB.r.vec vector of number of trials
#' @param NB.prob.vec vector of probabilities
#'
#' @return vector of binary probability, dichotomous threshold
#' @examples
#' NB.r.vec <- c(10, 3, 16)
#' NB.prob.vec <- c(0.65, 0.4, 0.88)
#'
#' # Compute binary probabilities, PMFs, and thresholds
#' p <- calc.bin.prob.NB(NB.r.vec, NB.prob.vec)
#'
#' @export

calc.bin.prob.NB <- function(NB.r.vec, NB.prob.vec){
  validation.NBparameters(NB.r.vec, NB.prob.vec)
  J <- length(NB.r.vec)
  p <- numeric(J)
  prop <- list(mode = vector)
  Mlocation <- numeric(J)

  for(i in 1:J){
    j = rnbinom(100000, size = NB.r.vec[i], prob = NB.prob.vec[i])
    j.table = table(j)/100000
    nb.upper.limit = as.numeric(names(j.table)[length(j.table)])

    while(  dnbinom(nb.upper.limit, size = NB.r.vec[i], prob = NB.prob.vec[i]) >= 1e-10){
      nb.upper.limit = nb.upper.limit +1
    }

    X.prop = dnbinom(c(0:nb.upper.limit) , size = NB.r.vec[i], prob = NB.prob.vec[i])
    names(X.prop) = as.character (c(0:nb.upper.limit))
    X_median = qnbinom(0.5, size = NB.r.vec[i], prob = NB.prob.vec[i])
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
