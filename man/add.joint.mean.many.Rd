\name{add.joint.mean.many}
\alias{add.joint.mean.many}
\alias{rem.joint.mean.many}
\alias{timeave.joint.mean.many}
\alias{add.joint.mean.one}
\alias{rem.joint.mean.one}
\alias{timeave.joint.mean.one}
\alias{add.cond.mean.many}
\alias{rem.cond.mean.many}
\alias{timeave.cond.mean.many}
\alias{add.cond.mean.one}
\alias{rem.cond.mean.one}
\alias{timeave.cond.mean.one}
\alias{hold.cond.mean.one}
\alias{add.joint.meanSq.one}
\alias{add.cond.meanSq.one}
\alias{addrem.joint.mean.one}
\alias{addrem.cond.mean.one}
\alias{addhold.joint.mean.one}
\alias{addhold.cond.mean.one}
\alias{remhold.joint.mean.one}
\alias{remhold.cond.mean.one}
\alias{add.joint.meanSq.one}
\alias{add.cond.meanSq.one}
\alias{rem.joint.meanSq.one}
\alias{rem.cond.meanSq.one}
\alias{hold.joint.meanSq.one}
\alias{hold.cond.meanSq.one}
\alias{all.cond.mean.PO}
\alias{all.cond.mean2.PO}
\title{Mean counts and particle time averages for birth-death processes
  with immigration}
\description{
  A set of functions for calculating the joint and conditional mean sufficient
  statistics for partially observed birth-death process with
  immigration. The sufficient statistcs are the number of
  births and immigrations, the mean number of deaths, and
  the time average of the number of particles.
  The conditional expectations
  of these quantities are calculated for a finite time interval,
  conditional on the number of particles at the beginning and the
  end of the interval.
}
\usage{
add.joint.mean.many(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024)
rem.joint.mean.many(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024)
timeave.joint.mean.many(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024)
add.cond.mean.many(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024,
 prec.tol=1e-12, prec.fail.stop=TRUE)
rem.cond.mean.many(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024,
prec.tol=1e-12, prec.fail.stop=TRUE)
timeave.cond.mean.many(t,lambda,mu,nu=0,X0=1,delta=0.001,n=1024,
prec.tol=1e-12, prec.fail.stop=TRUE)
add.joint.mean.one(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4)
rem.joint.mean.one(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4)
timeave.joint.mean.one(t,lambda,mu,nu=0,X0=1,Xt,delta=0.001,n=1024,r=4)
add.cond.mean.one(t,lambda,mu,nu=0,X0=1,Xt,trans.prob=NULL,
joint.mean=NULL,delta=1e-04,n=1024, r=4,
prec.tol=1e-12, prec.fail.stop=TRUE)
rem.cond.mean.one(t,lambda,mu,nu=0,X0=1,Xt,trans.prob=NULL,
joint.mean=NULL,delta=1e-04,n=1024, r=4,
prec.tol=1e-12, prec.fail.stop=TRUE)
timeave.cond.mean.one(t,lambda,mu,nu=0,X0=1,Xt,trans.prob=NULL,
joint.mean=NULL, delta=1e-04,n=1024,r=4,
prec.tol=1e-12, prec.fail.stop=TRUE)
hold.cond.mean.one(t,lambda,mu,nu=0,X0=1,Xt, trans.prob=NULL,joint.mean=NULL,
delta=1e-04,n=1024,r=4,prec.tol=1e-12, prec.fail.stop=TRUE)
add.joint.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt, joint.mean=NULL, delta = 0.001,
 n=1024,r=4)
add.cond.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt, trans.prob=NULL,
joint.mean=NULL, delta = 0.001, n
= 1024,r=4, prec.tol=1e-12, prec.fail.stop=TRUE)
addrem.joint.mean.one(t, lambda, mu, nu = 0, X0 = 1, Xt, delta = 0.001,
n = 1024,r=4)
addrem.cond.mean.one(t, lambda, mu, nu = 0, X0 = 1, Xt, trans.prob=NULL,
delta = 0.001,n = 1024, r=4,prec.tol=1e-12, prec.fail.stop=TRUE)
addhold.joint.mean.one(t, lambda, mu, nu = 0, X0 = 1, Xt, delta = 0.001,
n = 1024,r=4)
addhold.cond.mean.one(t, lambda, mu, nu = 0, X0 = 1, Xt, 
trans.prob=NULL, delta = 0.001, n = 1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE)
remhold.joint.mean.one(t, lambda, mu, nu = 0, X0 = 1, Xt, delta = 1e-04,
n = 1024,r=4)
remhold.cond.mean.one(t, lambda, mu, nu = 0, X0 = 1, Xt,
trans.prob=NULL, delta = 1e-04,
n = 1024,r=4, prec.tol=1e-12, prec.fail.stop=TRUE)
add.joint.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt, joint.mean=NULL,
delta = 0.001, n = 1024,r=4)
add.cond.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt,
trans.prob=NULL,joint.mean=NULL, delta = 0.001,
n= 1024, r=4, prec.tol=1e-12, prec.fail.stop=TRUE )
rem.joint.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt,
joint.mean=NULL, delta = 0.001, n = 1024,r=4)
rem.cond.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt, trans.prob=NULL,
joint.mean=NULL, delta = 0.001,n = 1024,r=4, prec.tol=1e-12, prec.fail.stop=TRUE)
hold.joint.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt, r=4,  n = 1024,
delta = 0.0001)
hold.cond.meanSq.one(t, lambda, mu, nu = 0, X0 = 1, Xt, trans.prob=NULL,
n= 1024,delta = 0.0001, r=4, prec.tol=1e-12, prec.fail.stop=TRUE)
all.cond.mean.PO(data,lambda,mu,nu=0,delta=0.001,n=1024, r=4,
prec.tol=1e-12, prec.fail.stop=TRUE)
all.cond.mean2.PO(data,lambda,mu,nu=0,delta=0.001,n=1024,r=4,
prec.tol=1e-12,  prec.fail.stop=TRUE)
}
\arguments{
  \item{t}{length of the time interval}
  \item{lambda}{per particle birth rate}
  \item{mu}{per particle death rate}
  \item{nu}{immigration rate}
  \item{X0}{starting state, a non-negative integer}
  \item{Xt}{ending state, a non-negative integer}
  \item{data}{CTMC_PO_1 or an analogous list. List isn't always
    accepted (in all.cond.mean functions it isn't). all.cond.means
    both accept CTMC_PO_1 or CTMC_PO_many.}
  \item{trans.prob}{Either NULL or a precomputed transition
    probability for a process with the parameters passed in.  This saves
    the repeated computation of the same transition probability
    for multiple conditional expectations over the same interval. If
    NULL, the probability will just be computed in the function.}
  \item{joint.mean}{This is a parameter in some of the computations for
    some squared means. It is either NULL or the corresponding (first-order,
    unsquared) mean.  If NULL the probability will just be
    computed in the function. (It is needed to convert a factorial mean
    to a squared mean.) Note that this is ALWAYS an unsquared mean,
    regardless of whether the function is a *.meanSq.* or a *.mean.*
    function.  In the latter case, if a non-NULL value is passed, the
    called function doesn't do much
    besides divide.}
  \item{delta}{increment length used in numerical differentiation}
  \item{n}{number of coefficients to pull off via FFT, in *.one
    functions this number determines the number of intervals in the
    Rieman sum approximation of the integral}
  \item{prec.tol}{"Precision tolerance"; to compute conditional means,
    first the joint means are computed and then they are normalized by
    transition probabilities.  The precision parameters govern the
    conditions under which the function will quit if these values are
    very small.  If the joint-mean is smaller than prec.tol then
    the value of prec.fail.stop decides whether to stop or continue.}
  \item{prec.fail.stop}{If true, then when joint-mean values are
    smaller than prec.tol the program stops; if false then it continues,
    usually printing a warning.}
  \item{r}{See numDeriv package; this is 'r' argument for
    grad/genD/hessian
    which determines how many richardson-method iterations are done.}
}
\details{
  
  Birth-death process is denoted by \eqn{X_t}
  
  Sufficient statistics are defined as

  \eqn{N_t^+} = number of additions (births and immigrations)
  
  \eqn{N_t^-} = number of deaths
  
  \eqn{R_t} = time average of the number of particles, \eqn{\int_0^t X_y dy}

  Function add.joint.mean.many returns a vector of length n, where the
  j-th element of the vector is equal to
  \deqn{E(N_t^+ 1_{X_t=j} | X_0=X0)}
    
  Function rem.joint.mean.many returns a vector of length n, where the
  j-th element of the vector is equal to
  \deqn{E(N_t^- 1_{X_t=j} | X_0=X0)}

  Function timeave.joint.mean.many returns a vector of length n, where the
  j-th element of the vector is equal to
  \deqn{E(R_t 1_{X_t=j} | X_0=X0)}

  Function add.cond.mean.many returns a vector of length n, where the
  j-th element of the vector is equal to
  \deqn{E(N_t^+ | X_0=X0, X_t = j)}
    
  Function rem.cond.mean.many returns a vector of length n, where the
  j-th element of the vector is equal to
  \deqn{E(N_t^- | X_0=X0, X_t=j)}

  Function timeave.cond.mean.many returns a vector of length n, where the
  j-th element of the vector is equal to
  \deqn{E(R_t | X_0=X0, X_t=j)}
  
  Function add.joint.mean.one returns
  \eqn{E(N_t^+ 1_{X_t=Xt} | X_0=X0)}

  Function rem.joint.mean.one returns
  \eqn{E(N_t^- 1_{X_t=Xt} | X_0=X0)}
  
  Function timeave.joint.mean.one returns
  \eqn{E(R_t 1_{X_t=Xt} | X_0=X0)}

  Function add.cond.mean.one returns
  \eqn{E(N_t^ | X_0=X0, X_t=Xt)}
  
  Function rem.cond.mean.one returns
  \eqn{E(N_t^- | X_0=X0, X_t=Xt)}
  
  Function timeave.cond.mean.one returns
  \eqn{E(R_t | X_0=X0, X_t=Xt)}

  Function add.joint.meanSq.one returns
  \eqn{E((N_t^-)^2, X_t=Xt | X_0=X0)}

  Function add.cond.meanSq.one returns
  \eqn{E((N_t^-)^2| X_0=X0, X_t=Xt )}
  
  Function addrem.joint.mean.one returns
  \eqn{E((N_t^- N_t^-) , X_t=Xt | X_0=X0)}
  
  Function addrem.cond.mean.one returns
  \eqn{E((N_t^- N_t^-)| X_0=X0, X_t=Xt )}
  
  all.cond.mean.PO and all.cond.mean2.PO compute the first
  and second order means respectively for a partially observed process
  (with possibly more than one observation point).  So they amalgamate
  the above functions and also apply them to multiple observations.
  The outcomes are labeled appropriately.

  Note that all.cond.mean.PO are not methods, they
  can accept either CTMC_PO_many or CTMC_PO_1 (via their use of
  CTMCPO2indepIntervals function).
  
  "Hold" and "timeave" are the same.
  
  The .many functions are less safe about differentiation right now.
  This should be changed in the future.
}
  
  
\author{Marc A. Suchard, Vladimir N. Minin, Charles Doss}

\examples{
library(DOBAD)
my.lambda = 2
my.mu = 3
my.nu =1
my.time = 0.5
my.start = 10
my.end = 2
my.n = 1024

#Calculate the mean number of additions (births and immigrations)
#conditional on "my.start" particles at time 0 and "my.end" particles at time "my.time"
add.cond.mean.one(t=my.time,lambda=my.lambda,mu=my.mu,nu=my.nu,X0=my.start,Xt=my.end)

#Calculate a vector mean number of deaths joint with "my.end" particles at
# time "my.time" and conditional on "my.start" particles at time 0
DOBAD:::rem.joint.mean.one(t=my.time,lambda=my.lambda,mu=my.mu,nu=my.nu,X0=my.start,Xt=my.end)

#Calculate a vector mean particle time averages conditional on
# "my.start" particles at time 0 and 1 to "my.n" particles at time "my.time" 
# WARNING: conditional expectations for large values of |X_0-X_t| may be
# unreliable
timeave.cond.mean.many(t=my.time,lambda=my.lambda,mu=my.mu,nu=my.nu,X0=my.start,n=my.n)[1:20]
}

\seealso{\code{\link{add.generator}}}


\keyword{math}
