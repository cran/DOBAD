\name{CTMC.simulate.piecewise}
\alias{CTMC.simulate.piecewise}
\title{
  Simulate from piecewise constant/homogeneous CTMC
}
\description{
  Via the CTMC.simulate function.
}
\usage{
CTMC.simulate.piecewise(rate.fns, jumpLim.fns, T.times, init.state)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rate.fns}{
    a LIST of rate functions corresponding to jumpLim.fns and T.times.
    Length is number of homogeneous pieces, we'll call it M.    
  }
  \item{jumpLim.fns}{
    a LIST of 'jumpLim' functions of length M like the list rate.fns.
    See the documentation for CTMC.simulate for an explanation of what
    each is.
  }
  \item{T.times}{
    Of length M+1 so that there are M intervals corresponding to rate.fns.
}
  \item{init.state}{
    A starting state for the simulated chain.
}
}

\value{
  An object of type CTMC.
}

\author{
Charles Doss
}