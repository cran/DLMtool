\name{Observation-class}
\Rdversion{1.1}
\docType{class}
\alias{Observation-class}
%%\alias{initialize,Observation-method}

\title{Class \code{"Observation"}}
\description{
An operating model component that controls the observation model
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Observation", OM)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Object of class \code{"character"} ~~ }
    \item{\code{ageMcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Cobs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Cbiascv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{CAAobs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{CALobs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Iobs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Mcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Kcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t0cv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Linfcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{LFCcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{LFScv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{B0cv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{FMSYcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{FMSY_Mcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{BMSY_B0cv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{rcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{A50cv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Dbiascv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Dcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Btbiascv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Btcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Fcurbiascv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Fcurcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{hcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Icv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{maxagecv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{beta}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Reccv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Irefcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Brefcv}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Crefcv}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Observation")}: ... }
	 }
}
\author{
T. Carruthers
}
\note{
Its questionable whether the hyperstability/hyperdepletion should be categorised as an observation model characteristic as it is most often driven by fleet dynamics (and therefore should be in the fleet object). Oh well its here and you might want to make it hyperstable beta < 1 or hyperdeplete beta > 1, only. 
}
\examples{
showClass("Observation")
}
\keyword{classes}
