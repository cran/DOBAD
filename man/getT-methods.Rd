\name{getT-methods}
\docType{methods}
\alias{getT}
\alias{getT-methods}
\alias{getT,BDMC-method}
\alias{getT,BDMC_many-method}
\alias{getT,CTMC-method}
\alias{getT,CTMC_many-method}
\alias{getT,CTMC_PO_1-method}
\alias{getT,CTMC_PO_many-method}
\title{ ~~ Methods for Function getT in Package `DOBAD' ~~}
\description{
 ~~ Methods for function \code{getT} in Package `DOBAD' ~~
}
\section{Methods}{
\describe{

\item{\code{signature(object = "BDMC")}}{
  Same as for CTMC.
}

\item{\code{signature(object = "BDMC_many")}}{
  Sum of time for each of component BDMCs.
}

\item{\code{signature(object = "CTMC")}}{
  time the chain is observed for.  Ie difference in first time we see
  the state and the last time.
}

\item{\code{signature(object = "CTMC_many")}}{
  Sum of time for components
}

\item{\code{signature(object = "CTMC_PO_1")}}{
  Difference in time first observation and last.
}

\item{\code{signature(object = "CTMC_PO_many")}}{
  Sum of time for components
}
}}
\keyword{methods}
