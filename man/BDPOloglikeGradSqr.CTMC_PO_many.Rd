\name{BDPOloglikeGradSqr.CTMC_PO_many}
\alias{BDPOloglikeGradSqr.CTMC_PO_many}
\title{
  Gradient-Squared of PartialData likelihood
}
\description{
  In Louis' 82 formula for the information of  partially observed data,
  the last term is the gradient-squared of the partial data likelihood.
  It doesn't have to be calculated because it's 0 at the MLE,
  but it's coded here for debugging purposes.
}
\usage{
BDPOloglikeGradSqr.CTMC_PO_many(partialDat, L, m, beta, n.fft = 1024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{partialDat}{
    CTMC_PO_many.
}
  \item{L}{
    lambda at which to calculate information; usually MLE.
}
  \item{m}{
    mu at which to calculate information; usually MLE.
}
  \item{beta}{
    known constant defining nu via nu=beta*lambda.
}
  \item{n.fft}{
    deprecated unused.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
