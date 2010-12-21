\name{Nplus}
\alias{Nplus}
\alias{Nminus}
\alias{Nplus.CTMC_PO_many}
\alias{Nminus.CTMC_PO_many}
\alias{holdTime}
\title{
  Calculate summary statistics for BDMCs and CTMC_POs
}
\description{
  Nplus (number jumps up), Nminus (number jumps down), and holdtime
  (waiting time weighted by the waiting state) are fundamental summary
  statistics for the Restricted Immigration BD model.  These functions
  compute those for BDMC, or compute the observed numbers for CTMC_PO_1
  or CTMC_PO_many.
}
\usage{
Nplus(sim)
Nminus(sim)
Nplus.CTMC_PO_many(ctmcpomany)
 Nminus.CTMC_PO_many(ctmcpomany)
 }

 \arguments{
   \item{sim}{
     Arg for Nplus, Nminus. BDMC generally.  Needs to have a getStates method.
   }
   \item{ctmcpomany}{A CTMC_PO_many.
   }
 }
 \details{
 }
 \value{
   Returns an integer, the number of jumps up.
 }
 \author{
   Charles Doss
 }


\examples{

Nplus(birth.death.simulant(1))
}
