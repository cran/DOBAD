\name{getBDjTimes}
\alias{getBDjTimes}
\title{
  Get Jump times of a BD process.
}
\description{
  get times of jumps, split into jumps up and jumps down.
}
\usage{
getBDjTimes(bdMC, getTimes = TRUE)
}
\arguments{
  \item{bdMC}{
    A BDMC
  }
  \item{getTimes}{
    Bool. If true returns times, otherwise returns indices of times vector.
}
}
\details{
  List with 2 components of times/indices.  First is times of jumps up,
  second is times of jumps down.
}
\value{
  If getTimes is TRUE:
  \item{timesup}{times of jumps up}
  \item{timesdown}{times of jumps down}

  If getTImes is FALSE:
  \item{indsup}{indices of times for jumps up}
  \item{indsdown}{indices of times for jumps down}
}

\author{
  charles doss
}
