\name{Stock-class}
\Rdversion{1.1}
\docType{class}
\alias{Stock-class}
%%\alias{initialize,Stock-method}

\title{Class \code{"Stock"}}
\description{
An operating model component that specifies the parameters of the population dynamics model
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Stock", OM)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Object of class \code{"character"} ~~ }
    \item{\code{maxage}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{R0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{M}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Msd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Mgrad}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{h}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Linf}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{K}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{t0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Ksd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Kgrad}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Linfsd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Linfgrad}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{recgrad}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{a}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{b}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{D}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ageM}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{ageMsd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Perr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Size_area_1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Frac_area_1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Prob_staying}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{Source}:}{Object of class \code{"character"} ~~ }
    \item{\code{AC}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{SRrel}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Stock")}: ... }
	 }
}
\author{
T. Carruthers
}
\examples{
showClass("Stock")
}
\keyword{classes}
