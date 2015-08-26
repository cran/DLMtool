# DLM_input MPs

matagelim<-function(x,DLM_data){ # Age at maturity is knife-edge vulnerability
  dependencies="DLM_data@AM, DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  Vuln<-1/(1+exp((DLM_data@AM[x]-(1:DLM_data@MaxAge))/(DLM_data@AM[x]*DLM_data@AM[x]*0.05)))
  c(Allocate, Effort, Spatial, Vuln)
}
class(matagelim)<-"DLM_input"

MRreal<-function(x,DLM_data){ # A Marine reserve in area 1 with spatial reallocation of effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(0,1)
  Vuln<-rep(NA,DLM_data@MaxAge)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRreal)<-"DLM_input"

MRnoreal<-function(x,DLM_data){ # A Marine reserve in area 1 with no spatial reallocation of effort
  dependencies="DLM_data@MaxAge"
  Allocate<-0
  Effort<-1
  Spatial<-c(0,1)
  Vuln<-rep(NA,DLM_data@MaxAge)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRnoreal)<-"DLM_input"

curE<-function(x,DLM_data){ # current effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-1
  Spatial<-c(1,1)
  Vuln<-rep(NA,DLM_data@MaxAge)
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE)<-"DLM_input"

curE75<-function(x,DLM_data){ #75% current effort
  dependencies="DLM_data@MaxAge"
  Allocate<-1
  Effort<-0.75
  Spatial<-c(1,1)
  Vuln<-rep(NA,DLM_data@MaxAge)
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE75)<-"DLM_input"
