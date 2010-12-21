\name{combineCTMC}
\Rdversion{1.1}
\alias{combineCTMC}
\title{
   Combine several CTMCs into one CTMC 
}
\description{
  Pastes together several CTMCs into one.  It doesn't check that the
  rules of the CTMCs are held to.  
}
\usage{
combineCTMC(sims)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sims}{
    a list each of whose element is a CTMC; so sims[[i]] is a CTMC. sims[[i]]
    can be of class "CTMC" or a list.
}
}
\details{
  Note that each CTMC should include "0" as its first time.  And the
  last state of sims[[i]] and the first state of sims[[i+1]] should
  "match" in that the user should check they follow the rules of
  whatever the generating process is for the CTMC.
}
\value{
  Returns a list (not a CTMC object!) with states, times, and T.
}
