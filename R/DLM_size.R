# DLM size methods

matsizlim<-function(x,DLM){
  return(1/(1+exp((DLM@AM[x]-(1:DLM@MaxAge))/(DLM@AM[x]*DLM@AM[x]*0.05))))
}
class(matsizlim)<-"DLM size"
