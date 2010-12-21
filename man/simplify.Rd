\name{simplify}
\Rdversion{1.1}
\alias{simplify}
\title{
  Transform Lists to Vectors
}
\description{
  Takes objects which are lists but are conceptually vectors, and
  transforms
  them into vector objects.
}
\usage{
simplify(simpleList)
}
\arguments{
  \item{simpleList}{
    A list each of whose components is a (numeric) vector of length 1.
  }
}
\details{
  simpleList is a list each of whose components is a (numeric) vector of
  length1; simplify returns a vectorized form of this list.
}
\value{
  numeric vector whose length is the number of components of simpleList.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
  The base R unlist function probably makes this redundant.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  unlist
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
\keyword{ }% __ONLY ONE__ keyword per line
