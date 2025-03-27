#' This function generates multivariate GPD data
#'
#' @param no.rows number of data
#' @param binObj intermediate correlation matrix object
#'
#' @return generated GPD data with user's specifications
#'
#' @export

genGPD <-function (no.rows, binObj)
{
  #browser()
  inter_bin = generate.binary(no.rows, binObj$pvec, binObj$intermat)
  Mydata = BinToGPD(binObj$pvec, binObj$GPDprop, binObj$Mlocation,
                    inter_bin)
  return(Mydata)
}
