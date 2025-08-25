#' Collapse binomial data outcomes to binary variable
#'
#' This function computes the binary probability and identifies a threshold
#' to split discrete binomial outcomes. It summarizes each binomial distribution by
#' providing the probability mass function, the probability of exceeding a median-based threshold,
#' and the location of the threshold for the binary split.
#'
#' @importFrom stats dbinom qbinom
#' @param B.n.vec Numeric vector of trial counts for each variable
#' @param B.prob.vec Numeric vector of success probabilities for each variable
#'
#' @return A list containing the binary probability vector, the list of binomial probability mass functions for each variable, and the corresponding threshold indices
#'
#' @examples
#' B.n.vec <- c(3, 4, 5)
#' B.prob.vec <- c(0.5, 0.5, 0.5)
#' p <- calc.bin.prob.B(B.n.vec, B.prob.vec)
#'
#' @export


calc.bin.prob.B <- function(B.n.vec, B.prob.vec){

  validation.Bparameters(B.n.vec, B.prob.vec)
  J <- length(B.n.vec)
  p <- numeric(J)
  prop <- list(mode = vector)
  Mlocation <- numeric(J)

  for(i in 1:J){
    b.upper.limit = B.n.vec[i]
    X.prop = dbinom(c(0:b.upper.limit) , size = B.n.vec[i], prob = B.prob.vec[i])
    names(X.prop) = as.character (c(0:b.upper.limit))
    X_median = qbinom(0.5, size = B.n.vec[i], prob = B.prob.vec[i])
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
