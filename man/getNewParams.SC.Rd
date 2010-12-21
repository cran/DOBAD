\name{getNewParams.SC}
\Rdversion{1.1}
\alias{getNewParams.SC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Solve for new parameters in restricted model in EM algorithm.
}
\description{
  Basically one step of the EM algorithm.  Given old parameters
  and the data, get the new parameters.
}
\usage{
getNewParams.SC(theData,oldParams, beta.immig,  dr = 0.001, r=4,n.fft =
1024, prec.tol, prec.fail.stop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oldParams}{
    Parameters from previous iteration
}
  \item{beta.immig}{
    immigrationrate = beta.immig * birthrate
}
  \item{theData}{
    The discretely observed BDI process.  Of class CTMC_PO_many, CTMC_PO_1, list.
}
  \item{dr}{
    tuning parameter for differentiation
  }
  \item{r}{Parameter for differentiation; see numDeriv package
    documentation.} 
  \item{n.fft}{
%%     ~~Describe \code{n.fft} here~~
}
  \item{prec.tol}{"Precision tolerance"; to compute conditional means,
    first the joint means are computed and then they are normalized by
    transition probabilities.  The precision parameters govern the
    conditions under which the function will quit if these values are
    very small.  If the joint-mean is smaller than prec.tol then
    the value of prec.fail.stop decides whether to stop or continue.}
  \item{prec.fail.stop}{If true, then when joint-mean values are
    smaller than prec.tol the program stops; if false then it continues,
    usually printing a warning.}
}
