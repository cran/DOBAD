\name{CTMCPO2indepIntervals}
\alias{CTMCPO2indepIntervals.CTMC_PO_1}
\alias{CTMCPO2indepIntervals.CTMC_PO_many}
\alias{CTMCPO2indepIntervals}
\title{
  Converts CTMC_PO (either CTMC_PO_1 or CTMC_PO_many) to independent intervals.
}
\description{
  The markov property means that conditional on endpoints, each
  interval of a markov chain is independent of the others.  
  For this reason   computations are often done on intervals.
}
\usage{
\method{CTMCPO2indepIntervals}{CTMC_PO_1}(partialDat)
\method{CTMCPO2indepIntervals}{CTMC_PO_many}(partialDat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{partialDat}{
    CTMC_PO_1 or CTMC_PO_many
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  This function converts data into a nx3 matrix where
  the first column is the starting state, the second is the ending state
  and the third is the length of time the interval spanned.  No
  distinction is made between data from "separate" units or
  separate intervals from the same markov chain.
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
