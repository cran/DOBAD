\name{getInitParams}
\Rdversion{1.1}
\alias{getInitParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get multiple starting parameters for EM
}
\description{
  This is an ad-hoc function that has some hardwired rules for grabbing
  a few starting parameter values given one passed in one.  You pass in
  your summary data and it makes a basic guess at starting parameters
  and then flips them and scales them in various ways depending on what
  numInitParams is.  It returns
  between 1 and 6 different values (pairs) of starting parameters.
  Useful to get more than one parameter
  for if you are automating the EM in some way.  Otherwise it
  just gives you the smart starting guess.
}
\usage{
getInitParams(numInitParams=1, summary.PO, T, beta.immig, diffScale)
}
\arguments{
  \item{numInitParams}{
    How many parameters you want returned; between 1 and 6.
    Note: the paramaters after the 1st are fairly arbitrary.
}
  \item{summary.PO}{
    Summary data from partially observed process.  "Nplus", "Nminus",
    and "Holdtime" should be names in that order of number of observed
    jumsp up, jumps down, and Holding time.
}
  \item{T}{
    total time of chain.
}
  \item{beta.immig}{
    Scalar multiple of lambda that gives you the immigration rate, ie
    immigrate = beta.immig * birthrate.
}
  \item{diffScale}{
    Note that we don't have a solution in the case mu == lambda.  So
    if the two are close then numerical differentiation requires smaller
    values essentially.  So usually pass something like "100*dr" where
    dr is the value that's passed through the add.joint.mean.*, etc
    (called delta) for numeric differentiation.
}
}
