\name{CTMC.simulate}
\alias{CTMC.simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate from ("regular") CTMC 
}
\description{
  Only CTMCs that have finite number of states to jump directly to
  starting from all given starting states are allowed.
}
\usage{
CTMC.simulate(rate.fn, jumpLim.fn, T.time, init.state)
}
\arguments{
  \item{rate.fn}{
    Rate function from N^2 -> R.  
}
  \item{jumpLim.fn}{
    Takes a state (integeR) as argument and
    returns an integer-pair. The 1st entry is the minimum
    possible state that can be jumped to from the argument
    as starting point, and the second is the maximum.
    These must be finite.
}
  \item{T.time}{
    length of time to simulate for.
}
  \item{init.state}{
    Starting state of sim.
}
}
\details{
  Simulates from a CTMC whose states are the integers.
  This version requires that each state can only jump
  to finitely many other states.  This information is
  encapsulated in jumpLim.fn.  This isn't fundamental
  but makes things proceed faster.
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
  Charles Doss
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
