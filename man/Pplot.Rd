\name{Pplot}
\alias{Pplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A projection by projection plot of F/FMSY and B/BMSY 
}
\description{
A shorter version of the plot method for MSEs that just shows the projected trends in stock status and over exploitation
}
\usage{
Pplot(MSEobj)
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
\examples{
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]]) # get data
Pplot(MSEtest1)
}