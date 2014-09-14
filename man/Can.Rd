\name{Can}
\alias{Can}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
What data-limited methods can be applied to this DLM data object?
}
\description{
An diagnostic tool that looks up the slot requirements of each method and compares this to the data available to limit the analysis to methods that have the correct data, do not produce errors and run within a time limit. Time limit is the maximum time taken to carry out five reps (stochastic samples) of a given method and is in units of seconds. 
}
\usage{
Can(DLM, timelimit = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DLM}{
A data-limited methods data object (class DLM)
}
  \item{timelimit}{
The maximum time (seconds) taken for a method to undertake 10 reps (this filters out methods that are too slow)
}
}
\examples{
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]]) # get data
avail('DLM')
Can(Red_snapper,0.1)
}
