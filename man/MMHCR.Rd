\name{MMHCR}
\alias{MMHCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mark Maunder's Harvest Control Rule
}
\description{
A harvest control rule that makes incremental adjustments to quota recommendations based on the apparent trend in surplus production
}
\usage{
MMHCR(x, DLM, reps = 100, alp = c(0.8, 1.2), bet = c(0.8, 1.2))
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
\item{alp}{
Condition for modifying the quota (bounds on change in abundance)
}
  \item{bet}{
Limits for how much the quota can change among years
}
}
\details{
Note that this isn't exactly what Mark has previous suggested and is stochastic in this implementation. 
}
\value{
A numeric vector of quota recommendations
}
\references{
http://www.iattc.org/Meetings/Meetings2014/MAYSAC/PDFs/SAC-05-10b-Management-Strategy-Evaluation.pdf
}
\author{
T. Carruthers
}
