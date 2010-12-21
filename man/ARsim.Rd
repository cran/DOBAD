\name{ARsim}
\alias{ARsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Accept-Reject Simulation
}
\description{
  Generic Code for acceptance-rejection sampling.
}
\usage{
ARsim(margSimFn, acceptFn, N, keepTestInfo = FALSE)
}
\arguments{
  \item{margSimFn}{
    This is a simulation function.  It should take no arguments (or
    have default values that will be used).  It should output one
    simulation; we will refer to it as being "type X".
  }
  \item{acceptFn}{
    Should take "type X" as argument and output True or False.
    It should actually output a list (T/F, extraInformation).
    ExtraInformation is generally whatever is being used to do the
    accept/reject part; it is technically only required if
    "keepTestInfo" is passed a True.
  }
  \item{N}{
    How many simulations total to run (regardless of eventual
    acceptance ratio); There is not currently a parameter for
    choosing to stop after a given number of acceptances.
  }
  \item{keepTestInfo}{
    True or False; if True then the result will be two lists, the second
    of which has the second output from acceptFn; generally the data
    used to decide whether to accept or reject the simulations.
  }
}
\details{
  Does accept-reject simultion; margSimFn is run N times; acceptFn
  decides which to keep and which to remove; 
}
\value{
  Returns a list with one (if keepTestInfo==FALSE) or two (if
  keepTestInfo==TRUE) components. The first is $acceptSims, and the second
  is $testVals.  The component acceptSims are the simulated values that
  were accepted.  To do further analysis, testVals is the corresponding
  list of information used to evaluate.
  
  In future:
  Will have option to pass all simulations as output, and to
  accept simulations (but presumably with a different acceptFn) to allow
  for more reuse.

  As an implementation note: want to do "replicate" inside this function
  so as to regulate the types of output.
}

\author{
  Charles Doss
}
