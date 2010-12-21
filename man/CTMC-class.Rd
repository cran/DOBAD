\name{CTMC-class}
\Rdversion{1.1}
\docType{class}
\alias{CTMC-class}

\title{Class "CTMC"}
\description{Continuous time Markov Chain class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CTMC", ...)}.

}
\section{Slots}{
	 \describe{
    \item{\code{states}:}{numerics; usually integers. }
    \item{\code{times}:}{numerics; an _increasing_ sequence. }
    \item{\code{T}:}{final "observation" time of the chain, or time at
      which it is posited to exist.  }
  }
}
\section{Methods}{
  \describe{
    \item{getStates}{\code{signature(object = "CTMC")}: ... }
    \item{getT}{\code{signature(object = "CTMC")}: ... }
    \item{getTimes}{\code{signature(object = "CTMC")}: ... }
  }
}

\author{Charles Doss }


\seealso{
  \link{getT,CTMC-method}
}
\examples{
showClass("CTMC")
}
\keyword{classes}
