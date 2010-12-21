\name{BD.EMInference.prodExpecs}
\alias{getBDsummaryExpecs}
\alias{getBDsummaryProdExpecs}
\title{Functions That Calculate Product Expectations Needed for
  Inference via EM Algorithm}
\description{
  In order to calculate the information matrix for partial data, several
  conditional expectations of products of sufficient statistics are
  needed.  
}
\usage{
getBDsummaryExpecs(sims, fnc=function(x){x})
getBDsummaryProdExpecs(sims, getsd=FALSE)
}
\arguments{
  \item{sims}{A list of Birth-Death CTMCs.}
  \item{fnc}{A one argument function. It should be a function from Reals
    to Reals, capable of accepting a vector as its argument}
  \item{getsd}{ Also return estimate of standard deviations of the prods
  }
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).  We observe the process 
  at a finite  set of times over a time
  interval  [0,T].

  In order to calculate the information matrix for partial data, several
  conditional expectations of products of sufficient statistics are
  needed.  
  We have a method for  simulation conditional on the data,
  sim.condBD, which we use to
  estimate these.


  Generally for getting the information matrix after running the EM
  algorithm,  sim.condBD is called to simulate with the given parameters
  (estimates, usually), and the output sims are passed.  It is often
  important that the same set of sims are used to get all the results if
  the goal is to create an information matrix.
}

\value{
  getBDsummaryExpecs simply returns (an estimate of) E(fnc(Nt+)),
  E(fnc(Nt-)), and E(fnc(Rt)), where Nt+, Nt-, and Rt are the numbre of
  jumps up, the number of jumps down, and the total holding time,
  respectively.  They are returned in that order, also with labels
  "Nplus", "Nminus", and "Holdtime".
  
  getBDsummaryProdExpecs returns E(Nt+ * Nt-), E(Nt+ * Rt), and E(Nt- *
  Rt), in that order, also with the labels "NplusNminus",
  "NplusHoldtime", "NminusHoldtime".

  Returns another row of with corresponding
  standard deviations if getsd=TRUE.
}
\author{Charles Doss}

\seealso{\code{\link{add.joint.mean.many}},
  \code{\link{sim.condBD}}}
