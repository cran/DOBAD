\name{EMutilities}
\alias{M.step.SC}
\alias{E.step.SC}
\title{Functions related to implementing the EM algorithm on partially observed
  Birth-Death Chain}
\description{
  These are functions for the EM algorithm on a partially observed
  linear birth-death process where the immigration rate is a constant scalar
  times the birthrate.  The ".SC" suffix refers to this constraint ("SC"
  stands for "Special Case").
  
  E.step.SC performs the "Expectation step" and M.step.SC performs the
  maximization step.

  BDloglikelihood.PO computes the log likelihood of a partially observed
  birth-death process.
}
\usage{
M.step.SC(EMsuffStats, T,beta.immig)
E.step.SC(theData, oldParams, beta.immig,  dr=0.001, n.fft=1024,
                      r=4, prec.tol, prec.fail.stop)
}
\arguments{
  \item{EMsuffStats}{Vector with names "Nplus", "Nminus", and
    "Holdtime", which are the number of jumps up, number of jumps down,
    and the total holding time, respectively.
    These often come from the E.step.SC function.
  }
  \item{T}{total Time the chain was observed for (ie usually the last
    observation time).}
  \item{beta.immig}{Immigration rate is constrained to be a multiple of
    the birth rate.  immigrationrate = beta.immig * lambda where lambda
    is birth rate.
  }
  \item{oldParams}{Parameters with which to compute the expectation}
  \item{n.fft}{Number of terms to use in the fast fourier transform or
    the riemann integration  when
    using the generating functions to compute probabilities or
    joint expectations for the birth-death process.  See the
    add.joint.mean.many, etc, functions.}
  \item{theData}{Partially observed chain. Has components $times and $states
    where dat$states[i] is the state observed at time dat$times[i]. (No
    $T component needed).}
  \item{dr}{Parameter for numerical differentiation}
\item{r}{Parameter for differentiation; see numDeriv package
  documentation.} 
  \item{prec.tol}{"Precision tolerance"; to compute conditional means,
    first the joint means are computed and then they are normalized by
    transition probabilities.  The precision parameters govern the
    conditions under which the function will quit if these values are
    very small.  If the joint-mean is smaller than prec.tol then
    the value of prec.fail.stop decides whether to stop or continue.}
  \item{prec.fail.stop}{If true, then when joint-mean values are
    smaller than prec.tol the program stops; if false then it continues,
    usually printing a warning.
  }
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].

  E.step.SC computes the needed expectations for the EM algorithm.
  These are the expectations of the sufficient statistics, conditional
  on the data.  These expectations are computed with respect to the
  measure given by oldParams, i.e. the chain governed by oldParams.
  
  M.Step.SC maximizes the partial-data likelihood given the passed in
  expecatations of the sufficient statistics, to get the parameter
  iterates for the  next step of the EM algorithm.
  (This is easy when we are in the "Special Case"
  where immigration is constrained.)

  BDloglikelihood.PO computes the log likelihood of the passed in
  birth-death process.
}

\value{ M.step.SC returns a length 2 vector with first element
  lambda-hat and second element mu-hat, the respective maximizers of the
  likelihood.
  
  E.step.SC returns a vector with names "Nplus", "Nminus", and
  "Holdtime."

  BDloglikelihood.PO returns a real number, the log-likelihood of the data.
}
\author{Charles Doss}

\seealso{\code{\link{EM.BD.SC}}}

