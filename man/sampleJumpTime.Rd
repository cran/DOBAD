\name{sampleJumpTime2}
\Rdversion{1.1}
\alias{sampleJumpTime2}
\alias{p.i}
\alias{f.i}
\title{
  Functions for Simulating Conditionally the first Jump of a chain.
}
\description{
  Simulates the time of the first jump given that we know whether it's
  up or down and have observed the chain at some point.
}
\usage{
sampleJumpTime2(T, a, b, up = TRUE, L, m, nu)
p.i(T, a, b, up, L, m, nu, n.fft = 1024, subdivisions = 100)
f.i(t, T, a, b, up, L, m, nu, n.fft = 1024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
    Time of (First) observation; i.e. time at which we know the state of
    the chain.
  }
  \item{t}{ time between 0 and T at which to get the density of tau for
    f.i.
  }
  \item{a}{
    Starting state of the chain at time 0 (X_0=a).  a>=0.
}
  \item{b}{
    Given State of the chain at time T. (X_T=a).  b>=0.
}
  \item{up}{
    Boolean, telling whether the first jump is up (TRUE) or down (FALSE).
}
  \item{L}{
    Linear birth rate.
}
  \item{m}{
    Linear   death rate.
}
  \item{nu}{
    Immigration rate.
  }
  \item{subdivisions}{
    Parameter for numerical integration ("integrate" R function).
  }
  \item{n.fft}{
    Parameter for numerical riemann integration ("by hand").
  }
}
\details{
  Let tau be the time of the first jump (after time 0) and X_t is the
  chain at time t.
  
  Function sampleJumpTime2
  simulates the value of the first jump of a BDMC, conditional on
  some data.    
  What is given is the state of the BDMC at the beginning and end,
  where the end is time T, as well as whether
  the first jump is up or down.
  (To simulate the chain over the time from 0 to T, repeatedly
  call this function  alternatively with p.i)

  The Function p.i
  simulates whether the first jump is up or down, given the data.
  i.e. if up==true then this returns the probability
  [tau < T AND X_tau = a+1]
  and if up==false then it's
  [tau < T AND X_tau = a-1].

  The function f.i returns the "density" at t of tau, ie
  "P([tau ==t AND X_tau = a+1] | X0=a, Xt=b)"
  and if up==false then it's
  "P([tau ==t AND X_tau = a-1] | X0=a, Xt=b)".
  Note that it doesn't actually integrate to 1.  p.i(T) is the integral
  of f.i to time T.  f.i(.)/p.i(T) is actually a density on [0,T].
  If X_T != X_0 then we know the first jump is before time T.
  However, keep in mind the event of interest is that the first jump is
  up (down) and at time t; even if we know there will be a first jump
  down, that doesn't prove the first jump won't be up.
  In general, we have
  \eqn{\int^T_0 f.i(up) + \int^T_0 f.i(down) + P(first jump is after
    time T) = 1}.
  That is, \eqn{\int_0^t f.i(up)(s) ds} is the probability the first
  jump is before time t _and_ it is up (given that the chain starts at a
  and ends at b).
}
\value{
  A time (real number) between 0 and T.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Charles Doss
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{p.i}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
