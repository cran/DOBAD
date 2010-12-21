\name{getPartialData}
\Rdversion{1.1}
\alias{getPartialData}
\title{
  Get "partially Observed " Chain from a fully observed one.
}
\description{
  This effectively turns "Truth" into "data," ie it is passed
  a fully observed chain and returns only a partially observed one.
}
\usage{
getPartialData(observeTimes, CTMC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observeTimes}{
    Times at which CTMC is to "be observed" ie at which "data" is to be gathered.
}
  \item{CTMC}{
    A continuous time markov chain.
}
}
\details{
  Returns a CTMC_PO_1, ie discretely observed CTMC,
  from observing CTMC at observeTimes
}
\value{
  Returns CTMC_PO_1.
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
