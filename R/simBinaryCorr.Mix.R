#' Calculate Intermediate Binary Correlations for Mixed Variables
#'
#' This function computes intermediate binary correlations to generate a mixture of
#' Generalized Poisson (GPD), Negative Binomial (NB), and Binomial (B) variables with
#' a specified target correlation structure. It iteratively adjusts pairwise correlations
#' in binary space to approximate the desired correlation matrix for the mixed variables.
#'
#' @param GPD.theta.vec Numeric vector of theta parameters for GPD variables (or `NULL` if none).
#' @param GPD.lambda.vec Numeric vector of lambda parameters for GPD variables (must match length of `GPD.theta.vec`).
#' @param NB.r.vec Numeric vector of dispersion parameters (`r`) for NB variables (or `NULL` if none).
#' @param NB.prob.vec Numeric vector of success probabilities for NB variables (must match length of `NB.r.vec`).
#' @param B.n.vec Numeric vector of number of trials for Binomial variables (or `NULL` if none).
#' @param B.prob.vec Numeric vector of success probabilities for Binomial variables (must match length of `B.n.vec`).
#' @param CorrMat Correlation matrix (must be symmetric positive definite with dimensions matching total variables).
#' @param no.rows Integer specifying the number of rows (samples) to generate during intermediate binary sampling.
#' @param steps Numeric step size (default = 0.025) for correlation adjustment in later iterations.
#'
#' @return A list containing:
#'   \item{Mixprop}{List of proportions for each variable's binary components.}
#'   \item{intermat}{Intermediate correlation matrix for binary variables (adjusted to be positive definite if needed).}
#'   \item{Mlocation}{List of location parameters for each variable.}
#'   \item{pvec}{Vector of binary probabilities for each variable.}
#'
#' @details
#' The function first calculates binary probabilities and properties for each distribution family (GPD, NB, Binomial)
#' using helper functions `calc.bin.prob.GPD`, `calc.bin.prob.NB`, and `calc.bin.prob.B`. It then iteratively adjusts
#' pairwise correlations in binary space to match the target correlation structure, using a step size for convergence.
#' If the intermediate matrix is not positive definite, it is adjusted using `Matrix::nearPD`.
#'
#'
#' @examples
#' \dontrun{
#' # Example with 2 GPD, 1 NB, and 1 Binomial variable
#' GPD.theta <- c(1, 2)
#' GPD.lambda <- c(0.5, 0.3)
#' NB.r <- 10
#' NB.prob <- 0.2
#' B.n <- 5
#' B.prob <- 0.5
#' CorrMat <- matrix(c(1, 0.3, 0.2, 0.1,
#'                    0.3, 1, 0.4, 0.2,
#'                    0.2, 0.4, 1, 0.3,
#'                    0.1, 0.2, 0.3, 1), 4, 4)
#' result <- simBinaryCorr.Mix(
#'   GPD.theta.vec = GPD.theta,
#'   GPD.lambda.vec = GPD.lambda,
#'   NB.r.vec = NB.r,
#'   NB.prob.vec = NB.prob,
#'   B.n.vec = B.n,
#'   B.prob.vec = B.prob,
#'   CorrMat = CorrMat,
#'   no.rows = 1000
#' )
#' }
#'
#' @importFrom Matrix nearPD
#' @export

simBinaryCorr.Mix = function(GPD.theta.vec=NULL, GPD.lambda.vec=NULL, NB.r.vec=NULL, NB.prob.vec=NULL, B.n.vec=NULL, B.prob.vec=NULL,
                             CorrMat, no.rows, steps = 0.025){

  n.variable = length(c(GPD.lambda.vec, NB.r.vec, B.n.vec))
  if(dim(CorrMat)[2] != n.variable) stop("Incorrect dimension of specified correlation matrix")
  no.GPD.var = length(GPD.theta.vec)
  no.NB.var = length(NB.r.vec)


  # Compute binary probabilities for each family if not NULL
  GPD.p <- if (!is.null(GPD.theta.vec)) calc.bin.prob.GPD(GPD.theta.vec, GPD.lambda.vec) else list(p=NULL, prop=NULL, Mlocation=NULL)
  NB.p  <- if (!is.null(NB.r.vec))      calc.bin.prob.NB(NB.r.vec, NB.prob.vec)          else list(p=NULL, prop=NULL, Mlocation=NULL)
  B.p   <- if (!is.null(B.n.vec))       calc.bin.prob.B(B.n.vec, B.prob.vec)             else list(p=NULL, prop=NULL, Mlocation=NULL)


  # GPD.p = calc.bin.prob.GPD(GPD.theta.vec, GPD.lambda.vec)
  # NB.p = calc.bin.prob.NB(NB.r.vec, NB.prob.vec)
  # B.p = calc.bin.prob.B(B.n.vec, B.prob.vec)

  pvec = c(GPD.p$p, NB.p$p, B.p$p)
  prop = c(GPD.p$prop, NB.p$prop, B.p$prop)
  Mlocation = c(GPD.p$Mlocation, NB.p$Mlocation, B.p$Mlocation)
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
    prop.pair <- list(prop[[pair.temp[1]]], prop[[pair.temp[2]]])
    change = 1
    iteration = 0
    cat("calculating the intermediate binary correlations pair Sigma", pair.temp , "\n")
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
  print(l == TRUE)
  cat("\n required ", iteration.total, " iterations to calculate intermediate binary correlations. \n")
  return(list(Mixprop = prop, intermat = intermat, Mlocation = Mlocation, pvec = pvec))
}
