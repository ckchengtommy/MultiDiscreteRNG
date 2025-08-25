#' This function implements Step 2 of the algorithm
#' It calculates the intermediate binary correlations.
#'
#' @importFrom utils combn
#' @importFrom matrixcalc is.positive.definite
#' @param GPD.theta.vec vector of theta values
#' @param GPD.lambda.vec vector of lambda values
#' @param CorrMat specified Correlation matrix
#' @param no.rows number of observations for generating Multivariate Binary data
#' @param steps Fraction of difference between the current and target matrix to be added in each iteration.
#'
#' @return intermediate multivariate binary Correlation matrix
#' @export
simBinaryCorr.GPD<- function (GPD.theta.vec, GPD.lambda.vec, CorrMat, no.rows, steps = 0.025){
  p = calc.bin.prob.GPD(GPD.theta.vec, GPD.lambda.vec)
  pvec = p$p
  prop= p$prop
  Mlocation = p$Mlocation
  n.variable = length(GPD.theta.vec)
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
#    cat("calculating the intermediate binary correlations pair Sigma", pair.temp , "\n")
    while (sum(change > 0.001) > 0) {
      iteration = iteration + 1
      cat("iteration:", iteration, "\n")
      cat("delta:", sum(change), "\n")
      inter_bin = generate.binaryVar(no.rows, pvec.pair, del.next)
      Mydata = BinToGPD(pvec.pair, prop.pair, Mlocation.pair, inter_bin)

      if (iteration < 15) {
        del.next = del.next + (CorrMat.pair - Mydata$Corr) * 0.9
      }
      else {
        del.next = del.next + (CorrMat.pair - Mydata$Corr) * steps
      }
      #      print(del.next)
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
  print(l == TRUE)
  cat("\n required ", iteration.total, " iterations to calculate intermediate binary correlations. \n")
  return(list(GPDprop = prop, intermat = intermat, Mlocation = Mlocation, pvec = pvec))
}
