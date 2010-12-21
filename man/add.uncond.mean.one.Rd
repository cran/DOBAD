\name{add.uncond.mean.one}
\alias{add.uncond.mean.one}
\alias{rem.uncond.mean.one}
\alias{hold.uncond.mean.one}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ENplus, ENminus, Eholdtime, unconditional on ending state.
}
\description{
  ENplus, ENminus, Eholdtime, unconditional on ending state.  i.e.
  sum over j of Eij(Nplus), etc.  Expected number of total jumps
  up/down/holdtime over the given interval, conditional on starting state.
}
\usage{
add.uncond.mean.one(t, X0, lambda, mu, nu, delta = 0.001, r = 4)
rem.uncond.mean.one(t, X0, lambda, mu, nu, delta = 0.001, r = 4)
hold.uncond.mean.one(t, X0, lambda, mu, nu, delta = 0.001, r = 4)
}
\arguments{
  \item{t}{
    time
}
  \item{X0}{
    starting state
}
  \item{lambda}{
    birth rate
}
  \item{mu}{
    death rate
}
  \item{nu}{
    immigration rate
}
  \item{delta}{
    paramter for derivative.
}
  \item{r}{
    parameter for derivative.
}
}
\details{
  Uses generating functions.
}
\value{
  Each return a numeric.
}

\author{
  Charles Doss
}

