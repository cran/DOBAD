\name{sim.condBD}
\alias{sim.condBD}
\alias{sim.condBD.1}
\alias{sim.condBD.main}
\alias{sim.condBD.main.CTMC_PO_1}
\alias{sim.condBD.main.list}
\alias{sim.condBD.main.default}
\title{Simulate birth-death process, Conditionally
  upon observing its state at  finite, discrete
  time points}
\description{
  Functions for simulating a linear-birth-death
  process with birth parameter lambda, death
  parameter mu, and immigration parameter
  modelParams["n"]*lambda, conditional upon
  observing it at a finite discrete set of times over a finite time
  interval,
  [0,T].
}
\usage{
sim.condBD(bd.PO = list(states = c(5, 7, 3),times = c(0, 0.4, 1)), N = 1,
           L = 0.5, m = 0.7, nu = 0.4, n.fft = 1024, prevSims=NULL)
\method{sim.condBD.main}{CTMC_PO_1}(bd.PO=
      new("CTMC_PO_1", states=c(5,7,3), times=c(0,.4,1)),
      L=.5, m=.7, nu=.4,  n.fft=1024)
\method{sim.condBD.main}{list}(bd.PO=
      list(states=c(5,7,3), times=c(0,.4,1)),
      L=.5, m=.7, nu=.4,  n.fft=1024)
\method{sim.condBD.main}{default}(bd.PO=
      list(states=c(5,7,3), times=c(0,.4,1)),
      L=.5, m=.7, nu=.4,  n.fft=1024)
sim.condBD.1(T=1, a=5, b=7, L=.5, m=.7, nu=.4, n.fft=1024)
}
\arguments{
  \item{N}{Number of simulations/replications to do}
  \item{bd.PO}{_P_artially _O_bserved process, i.e., the data.  Needs
    to have "components" 'states', 'times', and 'T'.  Can be CTMC_PO_1 or a
    list.}
  \item{a}{Starting state of chain}
  \item{b}{Ending state of chain}
  \item{T}{Duration of chain}
  \item{L}{Lambda, linear birth rate parameter}
  \item{m}{mu, linear death rate parameter}
  \item{nu}{nu, immigration rate parameter.}
  \item{n.fft}{Number of terms to use in the fast fourier transform when
    using the generating functions to compute probabilities or
    joint expectations for the birth-death process.  See the
    add.joint.mean.many, etc, functions.}
  \item{prevSims}{A possibly-NULL list of previous simulation results
    which will be prepended to the current simulations. 
    }
}
\details{

sim.condBD, given discretely observed data from a chain,
simulates N birthdeath chains conditionally on the data.

sim.condBD.1 is the helper; it simulates one piece of the chain, given
a starting and ending state.  That process is repeated (via markov prop)
to simulate across many data points.  So it only takes arguments a,b,and
T rather than a "CTMC_PO_1" as its data.



  The method of simulating exactly is essentially that of Hobolth and
  Stone (2008).  Briefly: we can write out the density of the time of
  the first jump up or of the first jump down.  We can integrate it from 0
  to T  to
  compute the probability there is a jump up (down, respectively) in
  that time interval. Thus we can simulate whether or not there is a
  jump, and whether it is up or down.  Then using the above mentioned
  density, we can simulate the time at which it occurs.  For more details
  see Hobolth and Stone (2008).
}

\value{An object of class BDMC, ie a (linear) Birth-Death Markov Chain,
except for sim.condBD which returns a list of objects of class BDMC.}

\seealso{\code{\link{add.joint.mean.many}}}
\author{Charles Doss}
\source{
  Hobolth and Stone. (2008)
  Efficient Simulation from Finite-state, Continuous-Time Markov Chains
  with Incomplete Observations, submitted Annals of Applied Statistics.
}  