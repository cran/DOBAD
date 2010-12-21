\name{process.prob.one}
\Rdversion{1.1}
\alias{process.prob.one}
\alias{process.prob.many}
\title{
  Calculate transition probability for linear birth death process.
}
\description{
  Calculate transition probability for linear birth death process.
}
\usage{
process.prob.one(t, lambda, mu, nu = 0, X0 = 1, Xt,eps.t=1e-10,
eps.params=1e-10, n = -111)
process.prob.many(t, lambda, mu, nu = 0, X0 = 1, n = 1024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
    Time for transition.
}
  \item{lambda}{
    Linear birth rate 
}
  \item{mu}{
    linear death rate
}
  \item{nu}{
    immigration rate.
}
  \item{X0}{
    starting state.
}
  \item{Xt}{
    ending state.
  }
  \item{n}{
    Deprecated; for backwards compatibility.
  }
  \item{eps.t}{ One precision level below which the function
    switches to using the generating function
    instead of the Orthogonal Polynomial Solution to calculate
      transition probability.  Needed when the parameters or time
      are close to a boundary for which the OPS isn't defined.
    }
    \item{eps.params}{ Another precision level like eps.t.
    }
  }
\details{
  Calculates P(X_t=Xt |X_0=X0).
}
