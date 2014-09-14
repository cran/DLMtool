\name{area1MPA}
\alias{area1MPA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An marine reserve in area 1
}
\description{
A naive spatial control that changes the catches according to a vector e.g. c(0.4,1) reduces catch to 40 percent in area 1 and does not reduce catch in area 2. 
}
\usage{
area1MPA(x, DLM)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in data / simulation object DLM
}
  \item{DLM}{
A data limited methods data object
}
}
\author{
T. Carruthers
}
\note{
This is a spatial F-reduction method only. It cuts F in area 1 with no reallocation. 
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (x, DLM) 
return(c(0, 1)), class = "DLM space")
}
