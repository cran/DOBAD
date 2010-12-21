\name{list2CTMC}
\Rdversion{1.1}
\alias{list2CTMC}
\title{
  Convert a list representation of a CTMC to the class version
}
\description{
    Convert a list representation of a CTMC to the class version
}
\usage{
list2CTMC(aCTMC)
}
\arguments{
  \item{aCTMC}{
    A CTMC represented as a list.  Should have a "states", "times"
    vectors and a T numeric.
}
}

\value{
  Returns the same data but as an object of class CTMC.
}
