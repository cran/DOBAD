\name{getBDinform.PO}
\alias{getBDinform.PO}
\title{Get Information Matrix for MLE estimates on
  Partially Observed Linear Birth Death (_S_pecial _C_ase
  with constrained immigration)}
\description{
  Assume we have data that is the state at discrete
  time points of a linear birth-death process, which has immigration
  parameter constrained to be a known constant times the birth rate.  
  After using EM Algorithm for estimating rate parameters of a linear
  Birth-Death process, this function gives the information matrix
  associated.
}
\usage{
getBDinform.PO(partialData,Lhat,Mhat, beta.immig,delta=.001,
            n=1024,r=4, prec.tol=1e-12,prec.fail.stop=TRUE)
}
\arguments{
  \item{Lhat}{MLE for lambda, the birth rate.}
  \item{Mhat}{MLE for mu, the death rate.}
  \item{beta.immig}{Immigration rate is constrained to be a multiple of
    the birth rate.  immigrationrate = beta.immig * lambda where lambda
    is birth rate.}
  \item{partialData}{Partially observed chain. CTMC_PO_1 or CTMC_PO_many}
  \item{n}{n for riemann integral approximatoin.}
  \item{r, delta,prec.tol,prec.fail.stop}{see help for, say, all.cond.mean.PO}
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].  After running the EM algorithm to do estimation,
  this function returns the information to get, for instance, asymptotic
  CIs.

  See the Louis paper for the method.

  To calculate the information matrix, the expecatations of the
  products of the sufficient statistics, conditional on the data, are
  needed.  They are calculated by Monte Carlo, and N is the number of
  simulations to run.
}

\value{Symmetric 2x2 matrix; First row/column corresponds to lambda,
  second corresponds to mu}
\author{Charles Doss}


\source{
  Louis, T A. (1982). Finding the observed information matrix when using the EM
  algorithm. \emph{J. Roy. Statist. Soc. Ser. B}. 44 226-233.
  
}