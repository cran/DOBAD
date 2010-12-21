\name{getIthJumpTime}
\alias{getIthJumpTime}
\alias{getIthState}
\alias{getIthJumpTimes}
\title{
  Get the jump times from a CTMC.
}
\description{
  Get the time of the ith jump
}
\usage{
getIthJumpTime(CTMC, i)
getIthJumpTimes(timesList, i)
getIthState(CTMC,i)
}
\arguments{
  \item{CTMC}{ A CTMC.
  }
  \item{timesList}{List of positive numerics, each of which is the
    list of times from a CTMC.
  }
  \item{i}{ Positive integer. Which jump to get the time of.
    Need to know the CTMC(s) jumped at least i times!
  }
}
\details{
  Need to know the CTMC(s) jumped at least i times.
}
\value{
  getIthJumpTime returns a single positive numeric.
  getIthJumpTimes returns a vector of positive numerics.
  getIthState returns a nonnegative integer, the state.
}

\author{
  Charles Doss
}
