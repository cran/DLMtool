\name{DD4010}
\alias{DD4010}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delay - Difference Stock Assessment with UMSY and MSY leading coupled with a 40-10 rule
}
\description{
A simple delay-difference assessment that estimates the OFL using a time-series of catches and a relative abundance index. In this iteration a 40-10 rule is imposed over the OFL recommendation.  
}
\usage{
DD4010(x, DLM, reps = 100)
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
\value{
A numeric vector of quota recommendations
}
\references{
Method based on equations of Carl Walters
}
\author{
T. Carruthers
}
