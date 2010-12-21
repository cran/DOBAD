\name{CTMC.simulate.piecewise}
\alias{CTMC.simulate.piecewise}
\title{
  Simulate from piecewise constant/homogeneous CTMC
}
\description{
  Via the CTMC.simulate function.
}
\usage{
CTMC.simulate.piecewise(rate.fns, jumpLim.fns, T.times, init.state)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rate.fns}{
    a LIST of rate functions corresponding to jumpLim.fns and T.times.
    Length is number of homogeneous pieces, we'll call it M.    
  }
  \item{jumpLim.fns}{
    a LIST of 'jumpLim' functions of length M like the list rate.fns.
    See the documentation for CTMC.simulate for an explanation of what
    each is.
  }
  \item{T.times}{
    Of length M+1 so that there are M intervals corresponding to rate.fns.
}
  \item{init.state}{
    A starting state for the simulated chain.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of type CTMC.
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
