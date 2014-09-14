\name{DynF}
\alias{DynF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dynamic Fratio method
}
\description{
The Fratio method with a controller that changes the level of F according to the relationship between Surplus production and biomass. Ie lower F when dSP/dB is positive and higher F when dSP/dB is negative. 
}
\usage{
DynF(x, DLM, yrsmth=10, gg=2, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in a data-limited methods object
}
  \item{DLM}{
A data-limited methods object
}
   \item{yrsmth}{
The number of historical recent years used for smoothing catch and biomass data
}  
  \item{gg}{
A gain parameter that modifies F according to the gradient in surplus production with biomass
}
  \item{reps}{
The number samples of the quota
}
}
\details{
The method smoothes historical catches and biomass and then infers the relationship between surplus production and biomass (as suggested by Mark Maunder and Carl Walters). The approach then regulates a F based policy according to this gradient in which F may range between two different fractions of natural mortality rate. 

The core advantage is the quota(t) is not strongly determined by quota (t-1) and therefore errors are not propagated. The result is method that tends to perform alarmingly well and therefore requires debunking ASAP. 
}
\value{
A numeric vector of quota recommendations
}
\references{
Made-up for this package. 
}
\author{
T. Carruthers
}
