#' Compute intermediate binary correlations for multivariate negative binomial data
#'
#' This function implements Step 2 of the algorithm to calibrate the intermediate
#' latent-normal correlation matrix used to generate correlated binary
#' variables for negative binomial margins. For each pair of variables, it iteratively
#' updates the latent correlation so that, after (i) generating correlated binary data
#' via \code{generate.binaryVar} and (ii) mapping back to negative binomial outcomes via
#' \code{BinToNB}, the empirical correlation of the resulting NB pair matches the
#' user-specified target correlation in \code{CorrMat}. The calibrated pairwise latent
#' correlations are then assembled into a full intermediate matrix, which is adjusted
#' to be positive definite if needed.
#'
#' @importFrom utils combn
#' @importFrom matrixcalc is.positive.definite
#' @param NB.r.vec vector of number of trials
#' @param NB.prob.vec vector of probabilities
#' @param CorrMat specified correlation matrix
#' @param no.rows number of observations for generating multivariate binary data
#' @param steps fraction of difference between the current and target matrix to be added in each iteration.
#' @return intermediate multivariate binary Correlation matrix
#' @export
#' @examples
#' r.vec <- c(3, 5)
#' p.vec <- c(0.7, 0.5)
#'
#' M<- c(0.45, 0.45)
#' N <- diag(2)
#' N[lower.tri(N)] <- M
#' cmat<- N + t(N)
#' diag(cmat) <- 1
#'
#' # In real data simulation, no.rows should set to 100000 for accurate data generation
#' # in the intermediate step.
#' binObj = simBinaryCorr.NB(NB.r.vec = r.vec, NB.prob.vec = p.vec, CorrMat = cmat,
#' no.rows = 20000, steps= 0.025)

simBinaryCorr.NB<- function (NB.r.vec, NB.prob.vec, CorrMat, no.rows, steps = 0.025){
  p = calc.bin.prob.NB(NB.r.vec, NB.prob.vec)
  pvec = p$p
  prop= p$prop
  Mlocation = p$Mlocation
  n.variable = length(NB.r.vec)
  pair.combn = combn(n.variable, 2)
  n.corr = ncol(pair.combn)
  intermat.vec <- c()
  iteration.total = 0

  #####Compute the intermediate pair-wise correlation value from the specified correlation matrix
  for(i in 1:n.corr){
    pair.temp <- pair.combn[,i]
    CorrMat.pair <- diag(2)
    CorrMat.pair[2,1] <- CorrMat.pair[1,2] <- CorrMat[pair.temp[1], pair.temp[2]]
    del.next = CorrMat.pair
    pvec.pair <- c(pvec[pair.temp[1]], pvec[pair.temp[2]] )
    Mlocation.pair <- c(Mlocation[pair.temp[1]], Mlocation[pair.temp[2]])
    prop.pair <- list(prop[[pair.temp[1]]], prop [[pair.temp[2]]])
    change = 1
    iteration = 0
    message("calculating the intermediate binary correlations pair Sigma", pair.temp , "\n")
    while (sum(change > 0.001) > 0) {
      iteration = iteration + 1
      #message("iteration:", iteration, "\n")
      #message("delta:", sum(change), "\n")
      inter_bin = generate.binaryVar(no.rows, pvec.pair, del.next)
      Mydata = BinToNB(pvec.pair, prop.pair, Mlocation.pair, inter_bin)

      if (iteration < 20) {
        del.next = del.next + (CorrMat.pair - Mydata$Corr) * 0.9
      }
      else {
        del.next = del.next + (CorrMat.pair - Mydata$Corr) * steps
      }
      #print(del.next)
      change = abs(CorrMat.pair - Mydata$Corr)

    }
    iteration.total = iteration.total + iteration
    intermat.vec[i] <- del.next[1,2]
  }
  N <- diag(n.variable)
  N[lower.tri(N)] <- intermat.vec
  intermat <- N + t(N)
  diag(intermat) <- 1
  if (is.positive.definite(intermat) == FALSE) {
    #warning("Tetrachoric correlation matrix is not positive definite)")
    intermat = as.matrix(nearPD(intermat, corr = TRUE, conv.norm.type = "I", keepDiag = TRUE)$mat)
    intermat = (intermat + t(intermat))/2
  }
  l = is.positive.definite(del.next)
  #print(l == TRUE)
  message("\n required ", iteration.total, " iterations to calculate intermediate binary correlations. \n")
  return(list(NBprop = prop, intermat = intermat, Mlocation = Mlocation, pvec = pvec))
}
