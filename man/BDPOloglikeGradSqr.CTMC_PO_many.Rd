\name{BDPOloglikeGradSqr.CTMC_PO_many}
\alias{BDPOloglikeGradSqr.CTMC_PO_many}
\title{
  Gradient-Squared of PartialData likelihood
}
\description{
  In Louis' 82 formula for the information of  partially observed data,
  the last term is the gradient-squared of the partial data likelihood.
  It doesn't have to be calculated because it's 0 at the MLE,
  but it's coded here for debugging purposes.
}
\usage{
BDPOloglikeGradSqr.CTMC_PO_many(partialDat, L, m, beta, n.fft = 1024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{partialDat}{
    CTMC_PO_many.
}
  \item{L}{
    lambda at which to calculate information; usually MLE.
}
  \item{m}{
    mu at which to calculate information; usually MLE.
}
  \item{beta}{
    known constant defining nu via nu=beta*lambda.
}
  \item{n.fft}{
    deprecated unused.
}
}