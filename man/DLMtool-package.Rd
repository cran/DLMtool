\name{DLMtool-package}
\alias{DLMtool-package}
\alias{DLMtool}
\docType{package}
\title{
Data-Limited Methods Toolkit
}
\description{
Simulation testing and implementation of data-limited fishery stock assessment methods
}
\details{
\tabular{ll}{
Package: \tab DLMtool\cr
Type: \tab Package\cr
Version: \tab 1.34\cr
Date: \tab 2014-09-11\cr
License: \tab GPL-2\cr
Depends: \tab methods\cr
}
}
\author{
Tom Carruthers
<t.carruthers@fisheries.ubc.ca>
}
\references{
Carruthers, T.R., Punt, A.E., Walters, C.J., MacCall, A., McAllister, M.K., Dick, E.J., Cope, J. 2014. Evaluating methods for setting catch limits in data-limited fisheries. Fisheries Research. 153: 48-68. 
}
\keyword{data-limited}
\keyword{fisheries}
\keyword{stock assessment}
\keyword{data-poor}
\keyword{management strategy evaluation}
\examples{

# --- Application to real fishery data ---

mydata<-new('DLM')             # create a new DLM data object and define:
mydata@Year<-2001:2010         # years
mydata@Cat<-matrix((11:20)*10*runif(10,0.5,1.5),nrow=1) # make up some annual catches
mydata@Ind<-matrix(seq(1.1,0.9,length.out=10)*runif(10,0.5,1.5),nrow=1) 
mydata@Mort<-0.2               # instantaneous natural mortality rate
mydata@Abun<-1000              # current abundance estimate (biomass)
mydata@FMSY_M<-0.5             # guess of the ratio of FMSY to natural mortality rate
mydata@vbLinf<-200             # maximum length
mydata@vbK<-0.2                # von B growth coefficient k
mydata@LFC<-50                 # length at first capture
Can(mydata)                    # what methods can be applied to these data
Needed(mydata)                 # what data are needed to make the rest work
mydata<-getQuota(mydata)       # calculate quotas
plot(mydata)                   # plot them
mydata<-Sense(mydata,"Fratio") # conduct a sensitivity analysis for one of the methods


# --- Management strategy evaluation ---

for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]]) # get data
avail('Stock')                 # what stock objects are available?
avail('Fleet')                 # what fleet objects are available?
avail('Observation')           # what observation model objects are available?
OM<-new('OM',Rockfish, Generic_fleet, Imprecise_Biased) # define a new operating model
avail('DLM quota')                                      # available quota methods
plot(MSEtest1)                                          # plot the results
summary(MSEtest1)                                       # and summarize them in a table

}
