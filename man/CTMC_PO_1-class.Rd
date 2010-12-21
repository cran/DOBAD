\name{CTMC_PO_1-class}
\Rdversion{1.1}
\docType{class}
\alias{CTMC_PO_1-class}




\title{Class "CTMC_PO_1"}
\description{Partially observed CTMC.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CTMC_PO_1", ...)}.
Like CTMCs butdon't have an ending time; the final observation time
serves that purpose.  
}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{times}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{BDsummaryStats.PO}{\code{signature(dat = "CTMC_PO_1")}: ... }
    \item{getStates}{\code{signature(object = "CTMC_PO_1")}: ... }
    \item{getTimes}{\code{signature(object = "CTMC_PO_1")}: ... }
	 }
}

\author{Charles Doss}

\examples{
showClass("CTMC_PO_1")
}
\keyword{classes}
