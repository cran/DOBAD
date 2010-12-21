\name{graph.CTMC}
\Rdversion{1.1}
\alias{graph.CTMC}
\alias{graph.CTMC.PO}
\title{
  Plot  CTMCs (possibly Partially Observed)
}
\description{
  Plot in piecewise fashion the CTMCs.  If it is partially observed, it
  just plots it as if it were fully observed; i.e., the chain is
  pretended to continue in the same state until we see a jump.
}
\usage{
graph.CTMC(CTMC, filename = NA, height = 6, width = 4.5, xlab="time",
ylab="State", ...)
graph.CTMC.PO(CTMC, filename = NA, height = 6, width = 4.5,
type="l", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CTMC}{
    Either a fully observed CTMC or a partially observed one.
    Partially observed ones don't have a "T" and fully observed do.
}
  \item{filename}{
    filename string, or NA.
}
  \item{height}{
    Passed to trellis if filename isn't NA.
}
  \item{width}{
    Passed to trellis if filename isn't NA.
  }
  \item{xlab}{
    X label.
  }
  \item{ylab}{
    Y label
  }
  \item{type}{
    As in the plot parameter.
  }
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  If your data is S4 class, you can use the plot method.  
}
\value{

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
  See also the s4 methods written for the plot function \code{\link{plot-methods}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
