#' This function generates multivariate NB data
#'
#' @param no.rows number of data
#' @param binObj intermediate correlation matrix object
#'
#' @return generated NB data with user's specifications
#'
#' @export


genNB <-function (no.rows, binObj)
{
  inter_bin = generate.binaryVar(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToNB(binObj$pvec, binObj$NBprop, binObj$Mlocation,
                   inter_bin)
  return(Mydata)
}
