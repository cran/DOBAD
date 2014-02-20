\name{EM.BD.SC}
\alias{EM.BD.SC}
\alias{EM.BD.SC.1}
\title{Expectation-Maximization on Linear Birth Death (_S_pecial _C_ase
  with constrained immigration)}
\description{
    EM Algorithm for estimating rate parameters of a linear
  Birth-Death process, in which the data is the state at discrete
  time points
}
\usage{
EM.BD.SC(dat, initParamMat, tol = 1e-04, M = 30, beta.immig,  dr =
1e-07, n.fft = 1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE,
verbose=1, verbFile=NULL)
EM.BD.SC.1(dat,init.params, tol = 0.001, M = 30, beta.immig,  dr =
1e-07, n.fft = 1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE,
verbose=1, verbFile=NULL)
}
\arguments{
  \item{initParamMat}{\eqn{n}x2 matrix. Each row is an initial parameter
    setting. \eqn{n} is the number of times to run the full EM
    algorithm. On the \eqn{n}th time the initial "guess" of the
    lambda is initParamMat[n,1] and of mu it's initParamMat[n,2]. Used
    to automate starting at dispersed values to ensure global maximum.
    Frequently \eqn{n} is one.    
  }
  \item{init.params}{Vector of length two, first number is the first
    guess for lambda, second is the guess for mu.  This is like a single
    row from initParamMat.}
  \item{M}{Maximum number of iterations for (each) EM algorithm to run
    through. EM algorithm stops at Mth iteration.}
    \item{tol}{Tolerance for EM algorithm; when two iterations are
      within tol of each other the algorithm ends.  Algorithm also ends
      after M iterations have been reached.
      (note: One can debate whether 'tol' should refer to the estimates or
      to the actual likelihood.  here it is the estimates, though).
    }
  \item{beta.immig}{Immigration rate is constrained to be a multiple of
    the birth rate.  immigrationrate = beta.immig * lambda where lambda
    is birth rate.}
  \item{n.fft}{Number of terms to use in the fast fourier transform or
    the riemann integration 
    when
    using the generating functions to compute probabilities or
    joint expectations for the birth-death process.  See the
    add.cond.mean.many, etc, functions.}
  \item{dat}{Partially observed chain. Either of class "CTMC_PO_many"
    for several independent histories, of class "CTMC_PO_1" for one
    history, or a list with components $times and $states
    where dat$states[i] is the state observed at time dat$times[i]
    (ie, if it is a list then it is analogous to "CTMC_PO_1").}
  \item{dr}{Parameter for numerical differentiation}
  \item{r}{Parameter for numerical differentiation; see numDeriv package
    documentation.} 
  \item{prec.tol}{"Precision tolerance"; to compute conditional means,
    first the joint means are computed and then they are normalized by
    transition probabilities.  The precision parameters govern the
    conditions under which the function will quit if these values are
    very small.  If the joint-mean is smaller than prec.tol then
    the value of prec.fail.stop decides whether to stop or continue.}
  \item{prec.fail.stop}{If true, then when joint-mean values are
    smaller than prec.tol the program stops; if false then it continues,
    usually printing a warning.}
    \item{verbose}{ Chooses level of printing.  Increasing from 0, which
    is no printing.  
  }
  \item{verbFile}{ Character signifying the file to print to.
    If NULL just to standard output.
  }
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].  Runs EM algorithm to do maximum likelihood.

  EM.BD.SC will run the algorithm on multiple starting values and return
  the history for the best starting value.
  EM.BD.SC.1 only runs the algorithm for one starting value.  Otherwise
  they are the same.  
}

\value{Returns a \eqn{M+1}x2 matrix, the nth row being the estimators at
  the nth iteration.  The first
  column is for lambda (birth), the second for mu (death).  If tol is
  reached before M iterations then many of the rows will be empty, but
  the M+1st always contains the estimators.}
\author{Charles Doss}

\seealso{\code{\link{add.cond.mean.many}}}

