#' This function generates multivariate Mixed Discrete Data
#'
#' @param no.rows number of data
#' @param binObj intermediate correlation matrix object
#'
#' @return generated mixed discrete data with user's specifications
#'
#' @export


genMix = function (no.rows, binObj)
{
  #browser()
  inter_bin = generate.binaryVar(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToMix(binObj$pvec, binObj$Mixprop, binObj$Mlocation,
                    inter_bin)
  return(Mydata)
}
