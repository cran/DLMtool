\name{SPMSY}
\alias{SPMSY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Surplus Production MSY method
}
\description{
A method that uses Martell and Froese (2012) method for estimating MSY to determine the OFL. Since their approach estimates stock trajectories based on catches and a rule for intrinsic rate of increase it also returns depletion. Given their surplus production model predicts K, r and depletion it is straighforward to calculate the OFL based on the Schaefer productivity curve
}
\usage{
SPMSY(x, DLM, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A position in a data-limited methods data object
}
\item{DLM}{
A data-limited methods data object (class DLM)
}
  \item{reps}{
The number of samples of the quota taken for the calculation of the quota
}
}
\details{
Requires the assumption that catch is proportional to abundance. Occasionally the rule that limits r and K ranges does not allow r-K pairs to be found that lead to the depletion inferred by the catch trajectories. In this case this method widens the search. 
}

\references{
Martell, S. and Froese, R. 2012. A simple method for estimating MSY from catch and resilience. Fish and Fisheries. DOI: 10.1111/j.1467-2979.2012.00485.x
}
\author{
T. Carruthers
}
