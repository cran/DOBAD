\name{BDMC-class}
\Rdversion{1.1}
\docType{class}
\alias{BDMC-class}

\title{Class "BDMC"}
\description{	Birth-Death(-Immigration) CTMCs. Changes in state must
  be by 1 only. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BDMC", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{times}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{T}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CTMC}"}, directly.
}
\section{Methods}{
  \describe{
    \item{BDsummaryStats}{\code{signature(sim = "BDMC")}: ... }
    \item{getStates}{\code{signature(object = "BDMC")}: ... }
    \item{getTimes}{\code{signature(object = "BDMC")}: ... }
	 }
}
\references{}
\author{}
\note{  }

\seealso{
  \link{BDsummaryStats,BDMC-method},
  \link{getT,BDMC-method}
  
}
\examples{
showClass("BDMC")
}
\keyword{classes}
