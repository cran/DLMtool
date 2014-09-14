\name{DD}
\alias{DD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Delay - Difference Stock Assessment with UMSY and MSY leading
}
\description{
A simple delay-difference assessment that estimates the OFL using a time-series of catches and a relative abundance index. 
}
\usage{
DD(x, DLM, reps = 100)
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
Method based on equations of Carl Walters (bug him with questions and expect colourful responses) 
}
\author{
T. Carruthers
}
\note{
This DD model is observation error only. Similar to many other assessment models it depends on a whole host of dubious assumptions such as temporally constant productivity and proportionality between the abundance index and real abundance. These tend to be the biggest drivers of performance for this method. 
}
