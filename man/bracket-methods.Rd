\name{bracket-methods}%% NOTE cannot have a name of "[-methods" or it
%% won't be read
\docType{methods}
\alias{[,CTMC_many,ANY,ANY,ANY-method}
\alias{[,CTMC_PO_many,ANY,ANY,ANY-method}
\title{Methods for Function \code{[} in Package \pkg{DOBAD}}
\description{
  Methods for function \code{[} in package \pkg{DOBAD}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "CTMC_many", i = "ANY", j = "ANY", drop = "ANY")}}{
  Returns a CTMC_many object from the list of CTMCs indicated by the subscripts.
}

\item{\code{signature(x = "CTMC_PO_many", i = "ANY", j = "ANY", drop = "ANY")}}{
  Returns a CTMC_PO_many object from the list of CTMCs indicated by the subscripts.
}
}}
\keyword{methods}
