% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BinToB.R
\name{BinToB}
\alias{BinToB}
\title{This function implements Step 5 of the algorithm.
It converts the multivariate binary data back to the original Binomial scale}
\usage{
BinToB(prop.vec.bin, BProp, Mlocation, bin.data)
}
\arguments{
\item{prop.vec.bin}{vector of binary probabilities}

\item{BProp}{Binary proportion}

\item{Mlocation}{locations of median in the vector}

\item{bin.data}{generated multivariate binary data}
}
\value{
multivariate Binomial data and its correlation matrix
}
\description{
This function implements Step 5 of the algorithm.
It converts the multivariate binary data back to the original Binomial scale
}
