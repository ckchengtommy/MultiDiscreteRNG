#' Generates multivariate binomial data from binary parameters
#'
#'
#'
#' @param no.rows number of data
#' @param binObj intermediate correlation matrix object
#'
#' @return generated Binomial data with user's specifications
#'
#' @export

genB <-function (no.rows, binObj)
{
  inter_bin = generate.binaryVar(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToB(binObj$pvec, binObj$BProp, binObj$Mlocation,
                  inter_bin)
  return(Mydata)
}
