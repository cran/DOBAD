\name{BDMC_many-class}
\Rdversion{1.1}
\docType{class}
\alias{BDMC_many-class}

\title{Class \code{"BDMC_many"}}
\description{%
  A vector of BDMCs. Changes in the state for each element list should
  be by 1 only.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("BDMC_many", ...)}
  and supplying a "CTMC_many" object.
}
\section{Slots}{
  \describe{
    \item{\code{CTMCs}:}{Object of class \code{"list"}, a
      "CTMC_many" object}
  }
}
\section{Extends}{
Class "CTMC_many", directly.
}
\section{Methods}{
  No methods defined with class "BDMC_many" in the signature.
}
\examples{
showClass("BDMC_many")
}
\keyword{classes}
