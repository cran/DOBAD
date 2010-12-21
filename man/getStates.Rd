\name{getStates}
\alias{getStates}
\alias{getStates,CTMC-method}
\alias{getStates,BDMC-method}
\alias{getStates,CTMC_PO_1-method}
\title{
  Get list of jump states.
}
\description{
  Object accessor. 
}
\usage{
getStates(object)
}
\arguments{
  \item{object}{ A CTMC or generalization.  has a list of jump states.
}
}
\details{
  Gets list of states at each associated time.
}
\value{
  numeric vector, integer valued.
}

\author{
  Charles Doss
}
