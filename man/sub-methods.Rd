\name{sub-methods}
\docType{methods}
\alias{[,CTMC_many-method}
\alias{[,CTMC_PO_many-method}
\alias{[[,CTMC_many-method}
\alias{[[,CTMC_PO_many-method}
\title{  Subscripting CTMCs}
\description{
  Subscripting methods for CTMCs.
}
\section{Methods}{
\describe{


\item{\code{signature(x = "CTMC_many")}}{
  Gets x@CTMC[i]
}

\item{\code{signature(x = "CTMC_PO_many")}}{
  Gets x@BDMCsPO[i]
}


\item{\code{signature(x = "CTMC_many")}}{
  Gets x@CTMC[[i]]
}

\item{\code{signature(x = "CTMC_PO_many")}}{
    Gets x@BDMCsPO[[i]]
}
}}
\keyword{methods}
