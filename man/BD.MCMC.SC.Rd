\name{BD.MCMC.SC}
\alias{BD.MCMC.SC}
\title{MCMC on Linear Birth Death Process}
\description{
  Bayesian parameter estimation via Gibbs sampler MCMC on
  Linear Birth Death process, (_S_pecial _C_ase of constrained immigration)
  in which the data is the state at discrete
  time points.
}
\usage{
BD.MCMC.SC(Lguess, Mguess, beta.immig, alpha.L, beta.L, alpha.M, beta.M, 
    data, burnIn = 100, N = 1000, n.fft = 1024,
verbose=1, verbFile=NULL, simMethod=-1,...) 

}
\arguments{
  \item{Lguess}{Starting point for \eqn{\lambda}{lambda}}
  \item{Mguess}{Starting point for \eqn{\mu}{mu}}
  \item{beta.immig}{Immigration rate = beta.immig \eqn{* \lambda}{* lambda}.}
  \item{alpha.L}{Shape parameter for prior for \eqn{\lambda}{lambda}}
  \item{beta.L}{Rate parameter for prior for \eqn{\lambda}{lambda}}
  \item{alpha.M}{Shape parameter for prior for \eqn{\mu}{mu}}
  \item{beta.M}{Rate parameter for prior for \eqn{\mu}{mu}}
  \item{data}{Partially observed chain. Has components $times and $states
    where dat$states[i] is the state observed at time dat$times[i].}
  \item{N}{Number of iterations to run the MCMC for.}
  \item{burnIn}{Number of initial parameter estimates  to throw out.
    (So need burnIn << N.) Choose burnIn==0 throws nothing away.}
  \item{n.fft}{Number of terms to use in the fast fourier transform or
    the riemann integration     when
    using the generating functions to compute probabilities or
    joint expectations for the birth-death process.  See the
    add.joint.mean.many, etc, functions.}
  \item{verbose}{ Chooses level of printing.  Increasing from 0, which
    is no printing.  
  }
  \item{verbFile}{ Character signifying the file to print to.
    If NULL just to standard output.
  }
  \item{simMethod}{ Switch between using Accept-reject simulation and
    using the exact simulation method.  If -1, the function attempts to
    determine the best one of the two for the given parameters.
    Value of 0 fixes it at AR, and 1 fixes it at the exact method.
  }
  \item{...}{ Unused at this point.}
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].
  This runs MCMC to do parameter estimation.  The method is Gibbs
  sampling, by augmenting the state space to include the  the
  fully observed chain.  Then Gibbs sampling is performed
  using the the conditional simulation of  sim.condBD and the fact that,
  given the fully observed chain as data, independent gamma priors
  are conjugate priors, with independent posteriors.
}

\value{Returns a \eqn{N-burnIn}x2 matrix, the nth row being the
  estimators/samples at  the nth iteration.  The first
  column is for lambda (birth), the second for mu (death).}
\author{Charles Doss}

\seealso{\code{\link{add.joint.mean.many}}}

