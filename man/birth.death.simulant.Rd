\name{birth.death.simulant}
\alias{birth.death.simulant}
\title{Simulation of birth-death processes with immigration}
\description{
  A set of functions for simulating and summarizing birth-death simulations
}
\usage{
birth.death.simulant(t,X0=1,lambda=1,mu=2,nu=1, condCounts=NULL)
}
\arguments{
  \item{t}{length of the time interval}
  \item{lambda}{per particle birth rate}
  \item{mu}{per particle death rate}
  \item{nu}{immigration rate}
  \item{X0}{starting state, a non-negative integer}
  \item{condCounts}{is either null or a numeric vector with items named
    "Nplus" and "Nminus" (possibly from BDsummaryStats).  
  }
}
\details{
  Birth-death process is denoted by \eqn{X_t}

  Function birth.death.simulant returns a BDMC object.    
}
  
  
\author{Marc A. Suchard}

\examples{
my.lambda = 2
my.mu = 3
my.nu =1
my.time = 0.5
my.start = 10
my.end = 2
my.n = 2000

# simulate a birth death trajectory
my.simulant=birth.death.simulant(t=my.time,X0=my.start,lambda=my.lambda,mu=my.mu,nu=my.nu)
print(my.simulant)

# summarize the simulated trajectory
BDsummaryStats(my.simulant)

}

\seealso{\code{\link{add.joint.mean.many}},\code{\link{add.generator}}}


\keyword{math}
