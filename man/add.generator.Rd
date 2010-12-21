\name{add.generator}
\alias{add.generator}
\alias{rem.generator}
\alias{timeave.laplace}
\alias{hold.generator}
\alias{process.generator}
\alias{addrem.generator}
\alias{addhold.generator}
\alias{remhold.generator}
\alias{addremhold.generator}
\title{Generating functions for birth-death processes with immigration}
\description{
  A set of generating functions for sufficient
  statistics for partially observed birth-death process with
  immigration. The sufficient statistcs are the number of
  births and immigrations, the mean number of deaths, and
  the time average of the number of particles.
}
\usage{
add.generator(r,s,t,lambda,mu,nu,X0)
rem.generator(r,s,t,lambda,mu,nu,X0)
timeave.laplace(r,s,t,lambda,mu,nu,X0)
hold.generator(w,s,t,lambda,mu,nu,X0)
process.generator(s,time,lambda,mu,nu,X0)
addrem.generator(u, v, s, t, X0, lambda, mu, nu)
remhold.generator( v, w, s, t, X0, lambda, mu, nu)
addhold.generator( u, w, s, t, X0, lambda, mu, nu)
addremhold.generator( u, v, w, s, t, X0, lambda, mu, nu)
}
\arguments{
  \item{r,u,v,w}{dummy variable attaining values between 0 and 1.
    We use r for the single-argument generators and u,v,w for
    births,deaths, and holdtime for the multi-variable generators syntax,
    generally.}
  \item{s}{dummary variable attaining values between 0 and 1}
  \item{t,time}{length of the time interval}
  \item{lambda}{per particle birth rate}
  \item{mu}{per particle death rate}
  \item{nu}{immigration rate}
  \item{X0}{starting state, a non-negative integer}
}
\details{
  Birth-death process is denoted by \eqn{X_t}
  
  Sufficient statistics are defined as

  \eqn{N_t^+} = number of additions (births and immigrations)
  
  \eqn{N_t^-} = number of deaths
  
  \eqn{R_t} = time average of the number of particles, \deqn{\int_0^t X_y dy}
  
  Function add.generator calculates
  \deqn{H_i^+(r,s,t) = \sum_{n=0}^\infty \sum_{j=0}^\infty
    Pr(N_t^+=n,X_t=j | X_o=i) r^n s^j}

  Function rem.generator calculates
  \deqn{H_i^-(r,s,t) = \sum_{n=0}^\infty \sum_{j=0}^\infty
    Pr(N_t^-=n,X_t=j | X_o=i) r^n s^j}
  
  Function timeave.laplace calculates
  \deqn{H_i^*(r,s,t) = \sum_{j=0}^\infty \int_0^\infty e^{-rx}
    dPr(R_t \le x, X_t=j | X_o=i) s^j}

  Function processor.generator calculates
  \deqn{G_i(s,t) = \sum_{j=0}^\infty
    Pr(X_t=j | X_o=i) r^n s^j}

  Function addrem.generator calculates
  \deqn{H_i(u,v,s,t) = \sum_{j=0}^\infty \sum_{n_1=0}^\infty  \sum_{n_2=0}^\infty 
    Pr(X_t=j, N_t^+=n_1, N_t^-=n_2 | X_o=i) u^{n_1} v^{n_2} s^j}

  Function addhold.generator calculates
  \deqn{H_i(u,,w,s,t) = \sum_{j=0}^\infty \sum_{n1 \ge 0}
    u^n_1 \int_0^\infty e^{-rx}
    dPr(R_t \le x, N_t^+=n_1, X_t=j | X_o=i) s^j}

  Function remhold.generator is the same as addhold.generator
  but with N- instead of N+.
}
\value{
 Numeric value of the corresponding generating function.
}


\author{Marc A. Suchard, Charles Doss}


\seealso{\code{\link{add.joint.mean.many}}}

\keyword{math}
