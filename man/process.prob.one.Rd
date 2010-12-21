\name{process.prob.one}
\Rdversion{1.1}
\alias{process.prob.one}
\alias{process.prob.many}
\title{
  Calculate transition probability for linear birth death process.
}
\description{
  Calculate transition probability for linear birth death process.
}
\usage{
process.prob.one(t, lambda, mu, nu = 0, X0 = 1, Xt,eps.t=1e-10,
eps.params=1e-10, n = -111)
process.prob.many(t, lambda, mu, nu = 0, X0 = 1, n = 1024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
    Time for transition.
}
  \item{lambda}{
    Linear birth rate 
}
  \item{mu}{
    linear death rate
}
  \item{nu}{
    immigration rate.
}
  \item{X0}{
    starting state.
}
  \item{Xt}{
    ending state.
  }
  \item{n}{
    Deprecated; for backwards compatibility.
  }
  \item{eps.t}{ One precision level below which the function
    switches to using the generating function
    instead of the Orthogonal Polynomial Solution to calculate
      transition probability.  Needed when the parameters or time
      are close to a boundary for which the OPS isn't defined.
    }
    \item{eps.params}{ Another precision level like eps.t.
    }
  }
\details{
  Calculates P(X_t=Xt |X_0=X0).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
