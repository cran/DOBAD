\name{BDloglikelihood.PO}
\alias{BDloglikelihood.PO}
\alias{BDloglikelihood.PO.CTMC_PO_1}
\alias{BDloglikelihood.PO.CTMC_PO_many}
\alias{BDloglikelihood.PO.list}
\Rdversion{1.1}
\title{
  Calculate log likelihood of Partially Observed BD process
}
\description{
  Calculates the  log likelihood of a "partially observed
  birth-death-immigration process."  
  
}
\usage{
\method{BDloglikelihood.PO}{CTMC_PO_1}(partialDat, L, m, nu,  n.fft = 1024)
\method{BDloglikelihood.PO}{CTMC_PO_many}(partialDat, L, m, nu,  n.fft = 1024)
\method{BDloglikelihood.PO}{list}(partialDat, L, m, nu,  n.fft = 1024)
BDloglikelihood.PO(partialDat, L, m, nu, n.fft = 1024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
    lambda, birth rate.
}
  \item{m}{
    mu, death rate.
}
  \item{nu}{
    nu, Immigration rate.
}
  \item{partialDat}{
    Either of class "CTMC_PO_many", or of
    class "CTMC_PO_1" or the latter's analog in list form, ie
    a list with the two components "states" and "times" for
    the "list" and default versions of this method.
}
  \item{n.fft}{
    precision for riemann integration / fast fourier transform.
}
}
\details{
  Immigration can be arbitrary here.  Calculates likelihood of the
  b-d-i proces when it is observed at discrete timepoints.
}
\value{
  Real number.
}
\author{
charles doss
}

\examples{
library(DOBAD)
T=25;
L <- .3
mu <- .6
beta.immig <- 1.2;
initstate <- 17;

#generate process
dat <- birth.death.simulant(t=T, lambda=L, m=mu, nu=L*beta.immig, X0=initstate);
#"observe" process
delta <- 2
partialData <- getPartialData( seq(0,T,delta), dat);
#calculate the likelihood
BDloglikelihood.PO(partialDat=partialData, L=L, m=mu, nu=beta.immig*L);
}