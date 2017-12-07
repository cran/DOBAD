\name{EM.BD.SC.cov.1sv}
\alias{EM.BD.SC.cov.1sv}
\title{
  Expectation-Maximization on Linear Birth Death (and constrained
  Immigration) with Covariates
}
\description{
  EM algorithm for maximum likelihood estimation of rate parameters of
  linear Birth-Death-Immigration processes in which data is the  state
  at discrete time points, and one has covariates.
}
\usage{
EM.BD.SC.cov.1sv(BDMCs.PO, ZZ.LL, ZZ.MM, coefs.LL.init, coefs.MM.init,
tol = 1e-04, M = 30, beta.immig,
dr = 1e-07, n.fft = 1024, r = 4,
prec.tol = 1e-12, prec.fail.stop = TRUE,
nlmiterlim = 100, nlmgradtol = 1e-09, nlmstepmax = 0.5, verbose = 1, verbFile = NULL)
}

\arguments{
  \item{BDMCs.PO}{
    Data for doing estimation.  See dat argument to EM.BD.SC.  Of class
    "CTMC_PO_many" (several independent histories) or "CTMC_PO_1" (one history).
  }
  \item{ZZ.LL}{
    Covariates (design matrix) for predicting lambda. (Could be
    duplicates of ZZ.MM.) The model is:   log lambda = ZZ.LL \%*\%
    gamma.LL, for some coefficients gamma.LL.
}
  \item{ZZ.MM}{
    Covariates (design matrix) for  predicting mu.  (Could be duplicates
    of ZZ.LL.)
    The model is: log mu = ZZ.mu \%*\% gamma.mu, for some coefficients gamma.mu.
}
  \item{beta.immig}{
    Immigration rate is constrained to be a multiple of
    the birth rate.  immigrationrate = beta.immig * lambda where lambda
    is birth rate.
}
  \item{coefs.LL.init}{
    Initial linear coefficients determining lambda, the birth rate.
}
\item{coefs.MM.init}{
  Initial linear coefficients determining mu, the death rate.
}
\item{tol}{Tolerance for EM algorithm; when two iterations are
  within tol of each other the algorithm ends.  Algorithm also ends
  after M iterations have been reached.
  (note: One can debate whether 'tol' should refer to the estimates or
  to the actual likelihood.  here it is the estimates, though).
}
  \item{M}{
    Maximum number of iterations for (each) EM algorithm to run
    through. EM algorithm stops at Mth iteration.
}
\item{dr}{
  Parameter for numerical differentiation.
}
\item{n.fft}{
  Number of terms to use in the fast fourier transform or
  the riemann integration
  when
  using the generating functions to compute probabilities or
  joint expectations for the birth-death process.  See the
  add.cond.mean.many, etc, functions.
}
  \item{r}{
    Parameter for numerical differentiation; see numDeriv package
    documentation.
}
  \item{nlmiterlim}{
    For optim() call in M.step.SC.cov.
}
  \item{nlmgradtol}{
    For optim() call in M.step.SC.cov.
}
  \item{nlmstepmax}{
    For optim() call in M.step.SC.cov.
}
  \item{prec.tol}{
    "Precision tolerance"; to compute conditional means,
    first the joint means are computed and then they are normalized by
    transition probabilities.  The precision parameters govern the
    conditions under which the function will quit if these values are
    very small.  If the joint-mean is smaller than prec.tol then
    the value of prec.fail.stop decides whether to stop or continue.
}
  \item{prec.fail.stop}{
    If true, then when joint-mean values are
    smaller than prec.tol the program stops; if false then it continues,
    usually printing a warning.
}
  \item{verbose}{
    Chooses level of printing.  Increasing from 0, which
    is no printing.
}
  \item{verbFile}{
    Character signifying the file to print to.
    If NULL just to standard output.
}
}
\details{
  Assume we have a linear-birth-death
  process \eqn{X_t} with birth parameter \eqn{\lambda}{lambda}, death
  parameter \eqn{\mu}{mu}, and immigration parameter
  \eqn{\beta \lambda}{beta*lambda} (for some known, real
  \eqn{\beta}{beta}).
  We use a log linear model so that  log lambda = ZZ.LL \%*\%
    gamma.LL, for some coefficients gamma.LL and similarly
    log lambda = ZZ.MM \%*\%
    gamma.MM for some coefficients gamma.MM.

  We observe the process
  at a finite  set of times over a time
  interval  [0,T].  Runs EM algorithm to do maximum likelihood.
}
\value{
  Returns a list with elements coeffs.LL, an \eqn{pp.LL}x\eqn{M+1}
  matrix, and coeffs.MM, an \eqn{pp.MM}x\eqn{M+1} matrix where \eqn{M}
  is the number of EM iterations and pp.LL and pp.MM are the number of
  coefficients for Lambda and Mu respectively.  The M+1 columns gives the
  final estimators.
}
\author{
  Charles R. Doss.
}

\examples{

library(Matrix)
library(functional)



set.seed(1234)
mm <- 30; ## num individuals. arbitrary.
pp <- 2; ## num covariates, = HALF the number parameters
ZZ <- matrix(rnorm(mm*pp, -2, .5), nrow=mm, ncol=pp); ## arbitrary ...
ZZ.l1 <- apply(ZZ, 1, Compose(sum,abs))
coefs0.LL <- rnorm(pp, 0, 1)
ZZ <- (1/ZZ.l1)*ZZ ## will need |coefs.LL_j-coefs0.MM.j|< logKK / max( ||z_i||_1)
KK <- 2
diffs0 <- (rbeta(pp, 2,2)-1/2) * log(KK) ## want |lambda-mu| within a factor of KK
coefs0.MM <- coefs0.LL+diffs0;
coefs0 <- matrix(c(coefs0.LL, coefs0.MM), nrow=pp,ncol=2)
theta0 <- exp(ZZ \%*\% coefs0);
initstates <- rpois(mm, 3)+1
Ts <- abs(rnorm(mm,1,1)) / (theta0[,1]*initstates)
bb <- 1.1; ##beta
arg <- cbind(Ts,theta0, bb*coefs0.LL, initstates);
colnames(arg) <- NULL
BDMCs <- apply(arg, 1,
function(aa){birth.death.simulant(aa[1],aa[5], aa[2],aa[3],aa[4])})
t.obs <- apply(cbind(rpois(mm,2)+2, Ts), 1,
 function(aa){sort(runif(aa[1], 0, aa[2]))}) ##at least 2 observs
BDMCs.PO <- apply(cbind(t.obs,BDMCs), 1,
function(aa){getPartialData(aa[[1]],aa[[2]])})
BDMCs.PO <- new("CTMC_PO_many", BDMCsPO=BDMCs.PO);



#### Run the EM: (commented for speed for CRAN checks)
##emRes1 <- EM.BD.SC.cov.1sv(BDMCs.PO,
##                           ZZ.LL=ZZ, ZZ.MM=ZZ,
##                           coefs.LL.init=coefs0.LL, ##initialize at truth (which are not MLEs)
##                           coefs.MM.init=coefs0.MM,
##                           tol=1e-4,
##                           M=2, ## for speed; increase.
##                           beta.immig=bb,
##                           dr=1e-7, n.fft=1024, r=4,
##                          prec.tol=1e-12, prec.fail.stop=TRUE,
##                           verbose=1, verbFile="BD_EM_covariates_tutorial.txt")


}

