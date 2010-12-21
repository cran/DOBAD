\name{getSubMC}
\alias{getSubMC}
\title{
  Extract a Sub Markov Chain
}
\description{
  Create a new sub markov chain from a given one.
}
\usage{
getSubMC(CTMC, T)
}
\arguments{
  \item{CTMC}{
    A CTMC object
}
  \item{T}{
    Time to cut off  the given CTMC to form a new one.
}
}
\details{
  Creates a new CTMC identical to the given CTMC from time 0 to T.
}
\value{
  a CTMC.
}
