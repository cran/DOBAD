\name{num.deriv}
\Rdversion{1.1}
\alias{num.deriv}
\alias{genDoneSided}
\alias{hessianOneSided}
\title{
  Numerical Differentiation
}
\description{
  Numerical derivative of one-d function defined on R.
}
\usage{
num.deriv(ftn, var, delta = 0.001, ...)
genDoneSided(func, x, sides, method = "Richardson",
method.args = list(eps = 1e-04,
    d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07), r = 4,
    v = 2), ...)
hessianOneSided(func, x, sides, method = "Richardson", method.args = list(eps = 1e-04, 
     d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07), r = 4,                
     v = 2), ...)                                                                 

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ftn, func}{
    A (differentiable) function defined on R.
}
  \item{var,x}{
    Value(s) at which to differentiate. x can be vector.
  }
  \item{sides}{of length equal to x; +1 is differentiate from above
    -1 from below.}
  \item{method}{see numDeriv package docs}
  \item{method.args}{see numDeriv package docs}
  \item{delta}{
    Small number defining accuracy of numeric derivative.
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  See the 'numDeriv' package from whence the genD function and hessian
  function come.  The versions here are one-sided adapatations of the
  originals from that package.  
}
\value{
  Real number.
}
