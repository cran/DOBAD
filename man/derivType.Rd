\name{derivType}
\alias{derivType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Helper for getting means from generating functions
}
\description{
  Choose whether to do one-sided or two-sided differentiation.  The
  latter is more effective/less unstable but not always defined.
}
\usage{
derivType(L, mu, eps = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
}
\details{
  Getting the means of interest from generating functions involves
  differentiation which is usually done numerically.  The functions
  of interest are fully defined on one side of the point of interest but
  have limited (if any) definition on the other side of the point.  For
  instance, if lambda=mu then the generator for the process N+ is not
  defined for r>1.  If lambda and mu are close then the process is
  defined for r>1 but very close to 1.  The function derivType takes
  lambda and mu and an epsilon and decides whether that epsilon is small
  enough to do a two sided derivative with epsilon as "h" or if a one
  sided derivative is needed.
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
