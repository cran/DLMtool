\name{CompSRA}
\alias{CompSRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Age-composition-based estimate of current stock depletion given constant Z
}
\description{
Estimates an OFL based on a Stock Reduction analysis fitted to current age-composition data. 
}
\usage{
CompSRA(x, DLM, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods data object
}
  \item{DLM}{
A data-limited methods data object
}
  \item{reps}{
The number of stochastic samples of the quota recommendation
}
}
\author{
T. Carruthers
}
\note{
This was supposed to be a bit like Jim Thorson's method, which is likely to be published soonish eg late 2014 - early 2015 (a draft came past me at some point in August). 
}
