\name{power.coef.one}
\Rdversion{1.1}
\alias{power.coef.one}
\alias{power.coef.many}
\title{
  Gets coefficients of a power series..  
}
\description{
  Reads off coefficients of a power series.
}
\usage{
power.coef.one(power.series, n = 1000, k, ...)
power.coef.many(power.series, n = 1024, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power.series}{
    A function from C to C.  Note that its single argument must be named
    "s".  Should be a power series.
}
  \item{n}{
    Parameter for numerical riemann integration or for FFT.
}
  \item{k}{
    The coefficient to get.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  
}
\value{
  A Real number, the kth coefficient, for .one,
  or a vector of coefficients for .many.
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
\keyword{  }
\keyword{  }% __ONLY ONE__ keyword per line
