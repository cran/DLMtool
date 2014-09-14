\name{Kplot}
\alias{Kplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
KOBE plot: a projection by projection plot of F/FMSY and B/BMSY
}
\description{
A standard KOBE plot by each method that also shows the percentage of methods that ended up in each quadrant
}
\usage{
Kplot(MSEobj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MSEobj}{
An object of class MSE
}
}
\author{
T. Carruthers 
}
\note{
Apologies for the nauseating shading. 
}
\examples{
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]]) # get data
Kplot(MSEtest1)
}