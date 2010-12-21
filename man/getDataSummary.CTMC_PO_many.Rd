\name{getDataSummary.CTMC_PO_many}
\alias{getDataSummary.CTMC_PO_many}
\alias{getDataSummary}
\title{
  Calculate Some Summarizing Information for the Given Data
}
\description{
  Computes some summarizing statistics for a CTMC_PO_many object
  and returns them, possibly also saving them to a file.
}
\usage{
getDataSummary.CTMC_PO_many(dat, file = "dataSummary.rsav")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    Discretely Observed BDI process.
}
  \item{file}{
    Filename to save to. 
}
}
\details{
  See the function definition for the variable names used.  Saving and
  loading to/from a file seemed like the simplest approach.
}

\author{
charles doss
}
