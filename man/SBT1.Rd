\name{SBT1}
\alias{SBT1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SBT simple harvest control rule
}
\description{
A harvest control rule that makes incremental adjustments to quota recommendations based on the apparent trend in surplus production
}
\usage{
SBT1(x, DLM, reps = 100, yrsmth=10, k1=1.5, k2=3, gamma=1)
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
The number of samples of the quota
}
  \item{yrsmth}{
The number of years for evaluating trend in relative abundance indices
}
  \item{k1}{
Control parameter
}
  \item{k2}{
Control parameter
}
  \item{gamma}{
Control parameter
}
}
\details{
This isn't exactly the same as the proposed methods and is stochastic in this implementation. The method doesn't tend to work too well under many circumstances possibly due to the lack of 'tuning' that occurs in the real SBT assessment environment. You could try asking Rich Hillary at CSIRO about this approach. 
}
\references{
http://www.ccsbt.org/site/recent_assessment.php
}
\author{
T. Carruthers
}
