\name{BDsummaryStats}
\alias{BDsummaryStats.PO}
\alias{BDsummaryStats.PO,list-method}
\alias{BDsummaryStats.PO,CTMC_PO_many-method}
\alias{BDsummaryStats.PO,CTMC_PO_1-method}
\alias{NijBD}
\alias{NijBD.CTMC_many}
\alias{waitTimes}
\alias{BDsummaryStats,BDMC_many-method}
\alias{BDsummaryStats,BDMC-method}
\alias{BDsummaryStats,list-method}
\alias{BDsummaryStats}
\title{Get summary statistics for EM Algorithm on Linear Birth-Death Process}
\description{
  When passed in a birth-death markov chain, this extracts the summary
  statistics that are needed for computing the MLE (if immigration is a
  fixed known constant multiple of birth).

  That is, BDsummaryStats returns the counts of the total number of
  jumps up, the total number of jumps down, and the total holding/waiting time
  (\deqn{\sum_i d(i)*i}, where d(i) is time spent in state i).

  BDsummaryStats.PO does something similar, but for a partially observed
  process.

  NijBD takes a BD CTMC and calculates the number of jumps up and the
  number of jumps down.

  waitTimes takes a CTMC and calculates the waiting time.
}
\usage{
BDsummaryStats(sim)
BDsummaryStats.PO(dat)
NijBD(BDhist)
NijBD.CTMC_many(BDhists)
waitTimes(stateHist, timeHist, T)

}
\arguments{
  \item{sim}{A fully observed BDMC (or list with $states, $times, $T),
    or a BDMC_many.    }
  
  \item{dat}{Partially observed CTMC (list with $states, $times, $T), no
    "BD" restrictions on the structure of the chain.}
  
  \item{BDhist}{States of a BDMC; can be either a vector of states (each
    differing from its predecessor by 1) or a BDMC in list or class form.}
  \item{BDhists}{CTMC_many object}
  
  \item{stateHist}{Vector of states (integers). Corresponds to timeHist.}
  \item{timeHist}{Vector of times (reals).  Corresponds to stateHist, i.e.
    stateHist[i] is the state at and after timeHist[i].}
  \item{T}{Total time the chain was observed for.}
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].

  If the process is fully observed then to calculate the MLEs, we need
  the number of jumps up, down, and the total holding time.
  BDsummaryStats takes a
  BD CTMC and returns these three values, in a vector, with the names
  "Nplus" and "Nminus" for the number of jumps up and number of jumps
  down, respectively, and the name "Holdtime" for the total holding
  time.

  If the process is not fully observed, then these statistics aren't
  known.  (The EM algorithm is essentially trying to get a best-guess of
  these statistics).  BDsummaryStats.PO returns, rather, a very naive
  guess.  It pretends that the process is essentially fully observed and
  computes the statistics from that.  Note it's not the same as calling
  BDsummaryStats since a BD process has stipulations on its format that
  a partially observed BD process doesn't. The values are returned with
  the same naming convention as BDsummaryStats.  

  NijBD takes the list of states of a BD CTMC, and returns a 2x(n+1)
  matrix,  where  n is the maximum state the chain visits.
  NijBD(arg)[1,k] is the number of jumps down from state k-1, and
  NijBD(arg)[2,k] is the number of jumps up from state k-1.
  
  waitTimes takes any fully observed CTMC and returns a numeric vector
  of length n+1 where the maximum state passed in is n.  The ith entry
  is the waiting time in the i-1st state.  So
  seq(0, to=n, by=1) \%*\% waitTimes gives the total holding time.
  
}

\value{See details}
\author{Charles Doss}

\seealso{\link{BDMC-class}}

