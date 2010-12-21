\name{bdARsimCondEnd}
\alias{bdARsimCondEnd}
\alias{bdARsimCondEnd.1}
\title{
  Conditional Simulation of BD via Accept-Reject
}
\description{
  Simulates linear birth-death processes conditional on observing the
  end time (or a series of discrete observations), via simple accept
  reject (ie marginal simulation and accepting if it has the right end state).
}
\usage{
bdARsimCondEnd(Naccepted = NULL, Ntotal = NULL,Nmax=NULL, bd.PO = new("CTMC_PO_1", states = c(5, 7, 3), times = c(0, 0.4, 1)), L = 0.5, m = 0.7, nu = 0.4)
bdARsimCondEnd.1(Naccepted = NULL, Ntotal = NULL, Nmax=NULL, T = 1.02, L = 0.3, m = 0.4, nu = 0.1, a = 8, b = 9)
}
\arguments{
  \item{Naccepted}{
    Number of accepted sims to have at the end.
    Naccepted overrides Ntotal.  If you want to use Ntotal, Naccepted
    should be NULL.  Note that the number of sims will be >= Naccepted,
    probably not exactly equal to Naccepted.  
  }
  \item{Ntotal}{
    Number of marginal sims to do; no guarantee of how many sims you
    will get out, but a better guarantee of how long it will take.
    If it gets no sims, it returns list().
  }
  \item{Nmax}{
    Different than Ntotal; it works with Naccepted.  The function
    quits when either it has Naccepted sims or when it has
    done Nmax attempts.  If it hits the max,returns whatever has been
    simulated so far, possibly list() if nothing.
  }
  \item{T}{
    Length of time of the chain.
}
\item{L}{
    Linear Birth rate.
}
  \item{m}{
    Linear death rate.
}
  \item{nu}{
    Immigration rate.
  }
  \item{bd.PO}{
    For bdARsimCondEnd, this is a list of observations essentially;
    Either class "CTMC_PO_1" or the analogous list.
  }
  \item{a}{
    Starting state.
}
  \item{b}{
    Ending state (when you have just one observation).
}
}
\details{
  Outputs a list of BDMC objects.  If Naccepted is not NULL then the
  list will be at least Naccepted long. 
}
\value{
  List of BDMC objects.
}

\author{
Charles Doss
}

\examples{
bdARsimCondEnd.1(Naccepted=10);  #default parameters; simulates at least10.
bdARsimCondEnd.1(Ntotal=10);  #default parameters; maybe end with none.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
