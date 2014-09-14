\name{Fleet-class}
\Rdversion{1.1}
\docType{class}
\alias{Fleet-class}
%%\alias{initialize,Fleet-method}

\title{Class \code{"Fleet"}}
\description{
The component of the operating model that controls fishing dynamics
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Fleet", OM)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Object of class \code{"character"} ~~ }
    \item{\code{nyears}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Spat_targ}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{AFS}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{age05}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Vmaxage}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Fsd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Fgrad}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Fleet")}: ... }
	 }
}
\author{
T. Carruthers
}
\examples{
showClass("Fleet")
}
\keyword{classes}
