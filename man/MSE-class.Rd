\name{MSE-class}
\Rdversion{1.1}
\docType{class}
\alias{MSE-class}
%%\alias{initialize,MSE-method}
%%\alias{plot,MSE-method}

\title{Class \code{"MSE"}}
\description{
A Management Strategy Evaluation object that contains information about simulation conditions and performance of data-limited methods
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MSE", Name, nyears, proyears, nmeths, meths, nsim, OMtable, Obs, B_BMSYa, F_FMSYa, Ba, FMa, Ca, OFLa)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Object of class \code{"character"} ~~ }
    \item{\code{nyears}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{proyears}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nmeths}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{meths}:}{Object of class \code{"character"} ~~ }
    \item{\code{nsim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{OM}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{Obs}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{B_BMSY}:}{Object of class \code{"array"} ~~ }
    \item{\code{F_FMSY}:}{Object of class \code{"array"} ~~ }
    \item{\code{B}:}{Object of class \code{"array"} ~~ }
    \item{\code{FM}:}{Object of class \code{"array"} ~~ }
    \item{\code{C}:}{Object of class \code{"array"} ~~ }
    \item{\code{quota}:}{Object of class \code{"array"} ~~ }
    \item{\code{SSB_hist}:}{Object of class \code{"array"} ~~ }
    \item{\code{CB_hist}:}{Object of class \code{"array"} ~~ }
    \item{\code{FM_hist}:}{Object of class \code{"array"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "MSE")}: ... }
    \item{plot}{\code{signature(x = "MSE")}: ... }
	 }
}
\author{
T. Carruthers
}
\keyword{classes}
