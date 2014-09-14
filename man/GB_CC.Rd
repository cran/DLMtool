\name{GB_CC}
\alias{GB_CC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Geromont and Butterworth Constant Catch Harvest Control Rule
}
\description{
A harvest control rule that aims for MSY subject to imperfect information.
}
\usage{
GB_CC(x, DLM, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A position in data-limited methods data object
}
  \item{DLM}{
A data-limited methods data object
}
  \item{reps}{
The number of quota samples
}
}
\details{
Note that this is my interpretation of their approach and is now stochastic. Currently it is generalized and is not 'tuned' to more detailed assessment data which might explain why in some cases it leads to stock declines. 
}

\references{
Geromont, H.F. and Butterworth, D.S. 2014. Complex assessment or simple management procedures for efficient fisheries management: a comparative study. ICES J. Mar. Sci. doi:10.1093/icesjms/fsu017
}
\author{
T. Carruthers
}
