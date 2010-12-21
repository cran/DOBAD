\name{Nij}
\Rdversion{1.1}
\alias{Nij}
\title{
  Count transitions in a fully observed CTMC.
}
\description{
  Returns a matrix with counts of transitions
}
\usage{
Nij(CTMC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CTMC}{
    (Fully observed) CTMC.
}
}
\details{
  the (i,j) element is the number of transitions from state (i-1) to
  state (j-1) that were observed.
}
\value{
  numeric matrix(ncol=max(CTMC)+1, nrow=max(CTMC)+1) where max(CTMC)
  is the max state of the CTMC.
}
