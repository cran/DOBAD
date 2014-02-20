\name{getBDinform}
\alias{getBDinform.lost.SC.manual}
\alias{getBDinform.PO.SC.manual}
\alias{getBDinform.full.SC.manual}
\title{Helpers for Getting Information Matrix for MLE estimates on
  Partially Observed Linear Birth Death (_S_pecial _C_ase
  with constrained immigration)}
\description{
  Assume we have data that is the state at discrete
  time points of a linear birth-death process, which has immigration
  parameter constrained to be a known constant times the birth rate.  
  After using EM Algorithm for estimating rate parameters of a linear
  Birth-Death process, these functions compute matrices related to the
  information matrix.
}
\usage{
getBDinform.full.SC.manual(ENplus, ENminus, L, m)
getBDinform.lost.SC.manual(ENplus, ENminus, EHoldtime,
                             ENplusSq, ENminusSq, EHoldtimeSq,
                             ENplusNminus, ENplusHoldtime, ENminusHoldtime,
                             L, m, beta.immig, T)
getBDinform.PO.SC.manual(ENplus, ENminus, EHoldtime,
                             ENplusSq, ENminusSq, EHoldtimeSq,
                             ENplusNminus, ENplusHoldtime, ENminusHoldtime,
                             L, m, beta.immig, T)
}
\arguments{
  \item{L}{Lambda, birth rate}
  \item{m}{Mu, death rate}
  \item{beta.immig}{Immigration rate is constrained to be a multiple of
    the birth rate.  immigrationrate = beta.immig * lambda where lambda
    is birth rate.}
  \item{T}{Amount of time process is observed for; corresponds to time
    window over which all the expectations are computed.}
  \item{ENplus}{Expectation of the \eqn{N_T^+}{N^+}, the
    number of jumps up , conditional on    the data.}
  \item{ENminus}{Expectation of \eqn{N_T^-}{N^-},
    the number of jumps down, conditional on    the data.}
  \item{EHoldtime}{Expectation of \eqn{R_T^+}{R}, the total holdtime, conditional on
    the data.}
  \item{ENplusSq}{Expectation of \eqn{N_T^{+2}}{(N^+)^2},
    the square of the number of jumps up, conditional on    the data.}
  \item{ENminusSq}{Expectation of \eqn{N_T^{-2}}{(N^-)^2},
    the square of the number of jumps down, conditional on the data.}
  \item{EHoldtimeSq}{Expectation of \eqn{R_T^{2}}{(R)^2},
    the square of the total holdtime, conditional on the data.}
  \item{ENplusNminus}{Expectation of \eqn{N_T^+ N_T^-}{N^+ * N^-},
    the product of the number of jumps up and the number of jumps down, conditional on the data.}
  \item{ENplusHoldtime}{Expectation of \eqn{N_T^+ R_T}{N^+ * R},
    the product of the number of jumps up and the total holdtime, conditional on the data.}
  \item{ENminusHoldtime}{Expectation of \eqn{N_T^- R_T}{N^- * R},
    the product of the number of jumps down and the total holdtime, conditional on the data.}
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].  Can run the EM algorithm to do maximum likelihood.
  These functions are used to then compute pieces related to the
  information matrix.
  
  See equations 3.2 and 3.3 in the Louis paper for the notation.
  
  getBDinform.lost.SC.manual computes \eqn{I_{X|Y}}{I_x|y}.

  getBDinform.full.SC.manual computes \eqn{I_{X}}{I_x}.

  getBDinform.PO.SC.manual computes \eqn{I_{Y}}{I_y} (i.e. the difference
  between the other two functions).

  They have the "manual" suffix because the user passes in the
  expectations.  Some of them can be computed analytically by the
  methods in this package, but others cannot, so those are usually done
  by Monte Carlo (conditional on the data) simulation.

  NOTE: To make sure the answers are coherent, it is important to pass
  in expectations that are consistent with each other.  For instance,
  if the expectations ENplus, ENminus, and EHoldtime are computed
  analytically but  simulations are used to estimate the rest, then
  the results may be nonsense, because the values passed in were not
  necessarily feasible expectations all from the same measure.
}

\value{Symmetric 2x2 matrix; First row/column corresponds to lambda,
  second corresponds to mu}
\author{Charles Doss}

\source{
  Louis, T A. (1982). Finding the observed information matrix when using the EM
  algorithm. \emph{J. Roy. Statist. Soc. Ser. B}. 44 226-233.
}