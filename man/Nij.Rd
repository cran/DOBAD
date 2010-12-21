\name{Nij}
\Rdversion{1.1}
\alias{Nij}
\title{
  Count transitions in a fully observed CTMC.
}
\description{
  Returns a matrix with counts of transitions
}
\usage{
Nij(CTMC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CTMC}{
    (Fully observed) CTMC.
}
}
\details{
  the (i,j) element is the number of transitions from state (i-1) to
  state (j-1) that were observed.
}
\value{
  numeric matrix(ncol=max(CTMC)+1, nrow=max(CTMC)+1) where max(CTMC)
  is the max state of the CTMC.
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
