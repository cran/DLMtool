# DLM Source code
# Jan 2015
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
utils::globalVariables(c("R0","Mdb","mod","i","nareas","nsim","dFfinal"))

tiny<-1E-15

# generic class finder  e.g. available('data.frame') avalable('DLM')
getclass<-function(x,classy)inherits(get(x),classy)
#avail<-function(classy)ls(envir=package:DLMtool)[unlist(lapply(ls(envir=.GlobalEnv),getclass,classy=classy))]
avail<-function(classy){
  unique(c(ls('package:DLMtool')[unlist(lapply(ls('package:DLMtool'),getclass,classy=classy))],
    ls(envir=.GlobalEnv)[unlist(lapply(ls(envir=.GlobalEnv),getclass,classy=classy))]))
}

# Define classes --------------------------------------------------------------------------------------------------------------------------
setClass("DLM_data",representation(Name="character",Year="vector",Cat="matrix",Ind="matrix",Rec="matrix",t="vector",AvC="vector",
                                  Dt="vector",Mort="vector",FMSY_M="vector",BMSY_B0="vector",Cref="vector",
                                  Bref="vector",Iref="vector",
                                  AM="vector",LFC="vector",LFS="vector",CAA="array",
                                  Dep="vector",Abun="vector",vbK="vector",vbLinf="vector",
                                  vbt0="vector",wla="vector",wlb="vector",steep="vector",
                                  CV_Cat="vector", CV_Dt="vector",CV_AvC="vector",CV_Ind="vector",CV_Mort="vector",
                                  CV_FMSY_M="vector",CV_BMSY_B0="vector",CV_Cref="vector",CV_Bref="vector",
                                  CV_Iref="vector",CV_Rec="vector",
                                  CV_Dep="vector",CV_Abun="vector",
                                  CV_vbK="vector",CV_vbLinf="vector",CV_vbt0="vector",
                                  CV_AM="vector",CV_LFC="vector",CV_LFS="vector",    
                                  CV_wla="vector",CV_wlb="vector",CV_steep="vector",sigmaL="vector",
                                  MaxAge="vector",
                                  Units="character",Ref="numeric",Ref_type="character",
                                  Log="list",params="list",
                                  PosMPs="vector",MPs="vector",
                                  OM="data.frame",Obs="data.frame",
                                  TAC="array",TACbias="array",
                                  Sense="array",
                                  CAL_bins="numeric",
                                  CAL="array",MPrec="vector"))

setMethod("initialize", "DLM_data", function(.Object,stock="nada"){
  # run an error check here
  if(file.exists(stock)){
    dat <- read.csv(stock,header=F,colClasses="character") # read 1st sheet
    dname<-dat[,1]
    dat<-dat[,2:ncol(dat)]

    .Object@Name<-dat[match("Name", dname),1]
    .Object@Year <-as.numeric(dat[match("Year",dname),dat[match("Year",dname),]!=""])
    .Object@Cat <-matrix(as.numeric(dat[match("Catch",dname),dat[match("Catch",dname),]!=""]),nrow=1)
    .Object@Ind<-matrix(as.numeric(dat[match("Abundance index",dname),1:length(.Object@Year)]),nrow=1)
    .Object@Rec<-matrix(as.numeric(dat[match("Recruitment",dname),1:length(.Object@Year)]),nrow=1)
    .Object@t<-as.numeric(dat[match("Duration t",dname),1])
    .Object@AvC<-as.numeric(dat[match("Average catch over time t",dname),1])
    .Object@Dt<-as.numeric(dat[match("Depletion over time t",dname),1])
    .Object@Mort<-as.numeric(dat[match("M",dname),1])
    .Object@FMSY_M<-as.numeric(dat[match("FMSY/M",dname),1])
    .Object@BMSY_B0<-as.numeric(dat[match("BMSY/B0",dname),1])
    .Object@Cref<-as.numeric(dat[match("Cref",dname),1])
    .Object@Bref<-as.numeric(dat[match("Bref",dname),1])
    .Object@Iref<-as.numeric(dat[match("Iref",dname),1])
    .Object@AM<-as.numeric(dat[match("Age at 50% maturity",dname),1])
    .Object@LFC<-as.numeric(dat[match("Length at first capture",dname),1])
    .Object@LFS<-as.numeric(dat[match("Length at full selection",dname),1])
    .Object@Dep<-as.numeric(dat[match("Current stock depletion",dname),1])
    .Object@Abun<-as.numeric(dat[match("Current stock abundance",dname),1])
    .Object@vbK<-as.numeric(dat[match("Von Bertalanffy K parameter", dname),1])
    .Object@vbLinf<-as.numeric(dat[match("Von Bertalanffy Linf parameter", dname),1])
    .Object@vbt0<-as.numeric(dat[match("Von Bertalanffy t0 parameter", dname),1])
    .Object@wla<-as.numeric(dat[match("Length-weight parameter a", dname),1])
    .Object@wlb<-as.numeric(dat[match("Length-weight parameter b", dname),1])
    .Object@steep<-as.numeric(dat[match("Steepness", dname),1])
    .Object@sigmaL<-as.numeric(dat[match("Sigma length composition", dname),1])
  
    .Object@CV_Cat<-as.numeric(dat[match("CV Catch", dname),1])
    .Object@CV_Dt<-as.numeric(dat[match("CV Depletion over time t", dname),1])
    .Object@CV_AvC<-as.numeric(dat[match("CV Average catch over time t", dname),1])
    .Object@CV_Ind<-as.numeric(dat[match("CV Abundance index", dname),1])
    .Object@CV_Mort<-as.numeric(dat[match("CV M", dname),1])
    .Object@CV_Rec<-as.numeric(dat[match("CV Rec", dname),1])
    .Object@CV_FMSY_M<-as.numeric(dat[match("CV FMSY/M", dname),1])
    .Object@CV_BMSY_B0<-as.numeric(dat[match("CV BMSY/B0", dname),1])
    .Object@CV_Cref<-as.numeric(dat[match("CV Cref", dname),1])
    .Object@CV_Bref<-as.numeric(dat[match("CV Bref", dname),1])
    .Object@CV_Iref<-as.numeric(dat[match("CV Iref", dname),1])
    .Object@CV_Dep<-as.numeric(dat[match("CV current stock depletion", dname),1])
    .Object@CV_Abun<-as.numeric(dat[match("CV current stock abundance", dname),1])
    .Object@CV_vbK<-as.numeric(dat[match("CV von B. K parameter", dname),1])
    .Object@CV_vbLinf<-as.numeric(dat[match("CV von B. Linf parameter", dname),1])
    .Object@CV_vbt0<-as.numeric(dat[match("CV von B. t0 parameter", dname),1])
    .Object@CV_AM<-as.numeric(dat[match("CV Age at 50% maturity", dname),1])
    .Object@CV_LFC<-as.numeric(dat[match("CV Length at first capture", dname),1])
    .Object@CV_LFS<-as.numeric(dat[match("CV Length at full selection", dname),1])
    .Object@CV_wla<-as.numeric(dat[match("CV Length-weight parameter a", dname),1])
    .Object@CV_wlb<-as.numeric(dat[match("CV Length-weight parameter b", dname),1])
    .Object@CV_steep<-as.numeric(dat[match("CV Steepness", dname),1])
    .Object@MaxAge<-as.numeric(dat[match("Maximum age", dname),1])
    .Object@MPrec<-as.numeric(dat[match("MPrec", dname),1])
    
    if(length(grep("CAL",dname))>1){
      CAL_bins<-as.numeric(dat[match("CAL_bins",dname),dat[match("CAL_bins",dname),]!=""])
      nCAL<-length(CAL_bins)-1
      .Object@CAL_bins<-CAL_bins
      CALdat<-grep("CAL ",dname)
      .Object@CAL<-array(as.numeric(as.matrix(dat[CALdat,1:nCAL])),dim=c(1,length(CALdat),nCAL))
    }

    CAAy<-grep("CAA",dname)[1:length(grep("CAA",dname))]
    CAAa<-sum(dat[CAAy[1],]!="")
    if(!is.na(CAAa)){
      .Object@CAA<-array(as.numeric(as.matrix(dat[CAAy,1:CAAa])),dim=c(1,length(CAAy),CAAa))
    }

    .Object@Units<-dat[match("Units", dname),1]
    .Object@Ref<-as.numeric(dat[match("Reference TAC",dname),1])
    .Object@Ref_type<-dat[match("Reference TAC type",dname),1]
    .Object@Log[[1]]<-paste("Created:", Sys.time())
    .Object@params<-new('list')
    .Object@OM<-data.frame(NA)
    .Object@Obs<-data.frame(NA)
    .Object@TAC<-array(NA,dim=c(1,1,1))
    .Object@TACbias<-array(NA,dim=c(1,1,1))
    .Object@Sense<-array(NA,dim=c(1,1,1))
    .Object@PosMPs<-NA
    .Object@MPs<-NA

  }else{
    if(stock!="MSE"){
      if(!is.na(stock))print("Couldn't find specified csv file, blank DLM object created")
    }
  }

  # Default values -------------------------------------------------------------

  if(NAor0(.Object@CV_Cat)).Object@CV_Cat<-0.2
  if(NAor0(.Object@CV_Dt)).Object@CV_Dt<-0.25
  if(NAor0(.Object@CV_AvC)).Object@CV_AvC<-0.2
  if(NAor0(.Object@CV_Ind)).Object@CV_Ind<-0.2
  if(NAor0(.Object@CV_Mort)).Object@CV_Mort<-0.2
  if(NAor0(.Object@CV_FMSY_M)).Object@CV_FMSY_M<-0.2
  if(NAor0(.Object@CV_BMSY_B0)).Object@CV_BMSY_B0<-0.045
  if(NAor0(.Object@CV_Cref)).Object@CV_Cref<-0.2
  if(NAor0(.Object@CV_Bref)).Object@CV_Bref<-0.2
  if(NAor0(.Object@CV_Iref)).Object@CV_Iref<-0.2
  if(NAor0(.Object@CV_Rec)).Object@CV_Rec<-0.2
  if(NAor0(.Object@CV_Dep)).Object@CV_Dep<-0.25
  if(NAor0(.Object@CV_Abun)).Object@CV_Abun<-0.25
  if(NAor0(.Object@CV_vbK)).Object@CV_vbK<-0.1
  if(NAor0(.Object@CV_vbLinf)).Object@CV_vbLinf<-0.1
  if(NAor0(.Object@CV_vbt0)).Object@CV_vbt0<-0.1
  if(NAor0(.Object@CV_AM)).Object@CV_AM<-0.2
  if(NAor0(.Object@CV_LFC)).Object@CV_LFC<-0.2
  if(NAor0(.Object@CV_LFS)).Object@CV_LFS<-0.2
  if(NAor0(.Object@CV_wla)).Object@CV_wla<-0.1
  if(NAor0(.Object@CV_wlb)).Object@CV_wlb<-0.1
  if(NAor0(.Object@CV_steep)).Object@CV_steep<-0.2
  if(length(.Object@sigmaL)==0).Object@sigmaL<-0.2
  if(length(.Object@CAA)==0).Object@CAA<-array(NA,c(1,1,1))
  if(length(.Object@CAL)==0).Object@CAL<-array(NA,c(1,1,1))
  if(length(.Object@CAL_bins)==0).Object@CAL_bins<-1
  if(length(.Object@TAC)==0).Object@TAC<-array(1,c(1,1))
  if(length(.Object@TACbias)==0).Object@TACbias<-array(1,c(1,1))
  if(length(.Object@Sense)==0).Object@Sense<-array(1,c(1,1))
  
  .Object
})


NAor0<-function(x){
  if(length(x)==0)return(TRUE)
  if(length(x)>0)return(is.na(x[1]))
}  

cv<-function(x)  sd(x)/mean(x)
sdconv<-function(m,sd)(log(1+((sd^2)/(m^2))))^0.5        # get log normal standard deviation from transformed space mean and standard deviation
mconv<-function(m,sd)log(m)-0.5*log(1+((sd^2)/(m^2)))    # get log normal mean from transformed space mean and standard deviation
alphaconv<-function(m,sd)m*(((m*(1-m))/(sd^2))-1)
betaconv<-function(m,sd)(1-m)*(((m*(1-m))/(sd^2))-1)
trlnorm<-function(reps,mu,cv)return(rlnorm(reps,mconv(mu,mu*cv),sdconv(mu,mu*cv)))

DLMdiag<-function(DLM_data,command="available",reps=5,timelimit=1){
  funcs1<-c(avail("DLM_output"),avail("DLM_input"))
  good<-rep(TRUE,length(funcs1))
  report<-rep("Worked fine",length(funcs1))
  test<-new('list')
  timey<-new('list')
  options(show.error.messages = FALSE)
  for(y in 1:length(funcs1)){
    if(class(match.fun(funcs1[y]))=="DLM_output"){
      time1<-Sys.time()
      suppressWarnings({
        setTimeLimit(timelimit*1.5)
        test[[y]]<-try(do.call(funcs1[y],list(x=1,DLM_data=DLM_data,reps=5)),silent=T)
        setTimeLimit(Inf)
      })
    }else{
      time1<-Sys.time()
      suppressWarnings({
        setTimeLimit(timelimit*1.5)
        test[[y]]<-try(do.call(funcs1[y],list(x=1,DLM_data=DLM_data)),silent=T)
        setTimeLimit(Inf)
      })
    }
    time2<-Sys.time()
    timey[[y]]<-time2-time1
    if(class(test[[y]])=="try-error"){
      report[[y]]<-"Insufficient data"
      good[[y]]<-FALSE
    }else if(sum(is.na(test[[y]]))==length(test[[y]])){
      report[[y]]<-"Produced all NA scores"
      good[[y]]<-FALSE
    }
    if(timey[[y]]>timelimit){
      report[[y]]<-"Exceeded the user-specified time limit"
      good[[y]]<-FALSE
    }
  } # end of funcs
  options(show.error.messages = TRUE)
  if(command=="available")return(funcs1[good])
  if(command=="not available")return(cbind(funcs1[!good],report[!good]))
  if(command=="needed")return(needed(DLM_data,funcs=funcs1[!good]))
}

Required<-function(funcs=NA){
  if(is.na(funcs[1]))funcs<-c(avail("DLM_output"),avail("DLM_input"))
  slots<-slotNames('DLM_data')
  slotnams<-paste("DLM_data@",slotNames('DLM_data'),sep="")
  repp<-rep("",length(funcs))

  for(i in 1:length(funcs)){
    temp<-format(match.fun(funcs[i]))
    temp<-paste(temp[1:(length(temp))],collapse=" ")
    rec<-""
    for(j in 1:length(slotnams))if(grepl(slotnams[j],temp))rec<-c(rec,slots[j])
    if(length(rec)>1)repp[i]<-paste(rec[2:length(rec)],collapse=", ")
  }
  cbind(funcs,repp,deparse.level=0)
}

needed<-function(DLM_data,funcs=NA){
  if(is.na(funcs[1]))funcs<-avail("DLM_output")
  slots<-slotNames('DLM_data')
  slotnams<-paste("DLM_data@",slotNames(DLM_data),sep="")
  repp<-rep("",length(funcs))

  for(i in 1:length(funcs)){
    temp<-format(match.fun(funcs[i]))
    temp<-paste(temp[1:(length(temp))],collapse=" ")
    rec<-""
    for(j in 1:length(slotnams)){
      if(grepl(slotnams[j],temp)&NAor0(slot(DLM_data,slots[j])))rec<-c(rec,slots[j])
    }
    repp[i]<-paste(funcs[i],": ",paste(rec[2:length(rec)],collapse=", "),sep="")
  }
  repp
}

# Primary functions
Can<-function(DLM_data,timelimit=1)DLMdiag(DLM_data,"available",timelimit=timelimit)
Cant<-function(DLM_data,timelimit=1)DLMdiag(DLM_data,"not available",timelimit=timelimit)
Needed<-function(DLM_data,timelimit=1)DLMdiag(DLM_data,"needed",timelimit=timelimit)

demographic2=function(log.r,M,amat,sigma,K,Linf,to,hR,maxage,a,b){		#switch on and off to use either S or m in MC simulations
  r=exp(log.r)
  lx=exp(-M)^((1:maxage)-1)				#survivorship
  logNormDensity=(dnorm(x=log((1:maxage)),mean=log(amat),sd=sigma))/(1:maxage)	#Maturity ogive calculation
  logNormDensity[1]=0
  sumlogNormDen=sum(logNormDensity)
  NormalisedMaturity=logNormDensity/sumlogNormDen
  proportionMat[1]=NormalisedMaturity[1]
    for(i in 2:maxage)  proportionMat[i]=proportionMat[i-1]+NormalisedMaturity[i]
  TL=Linf*(1-exp(-K*((1:maxage)-to)))		#length at age
  Wa=a*TL^b					#wegith at age
  SurvWeiMat=lx*Wa*proportionMat	#survivorship X weight X maturity
  SBPR=sum(SurvWeiMat)		#Spawner biomass per recruit
  RPS=1/(SBPR*(1-hR)/(4*hR)) # Beverton Holt
  #RPS=(5*hR)^(5/4)/SBPR			# Ricker Recruitment per spawner biomass
  RPF=Wa*proportionMat*RPS		#Recruits per female
  Lotka=lx*RPF*exp(-(1:maxage)*r)
  sumLotka=sum(Lotka)
  epsilon=(1-sumLotka)^2	 				#objective function
  return(list(epsilon=epsilon,r=r))
}

demofn<-function(log.r,M,amat,sigma,K,Linf,to,hR,maxage,a,b)demographic2(log.r,M,amat,sigma,K,Linf,to,hR,maxage=maxage,a,b)$epsilon

proportionMat<-TL<-Wa<-SurvWeiMat<-r<-lx<-logNormDensity<-sumlogNormDen<-NULL
proportionMat=vector()


TAC<-function(DLM_data,MPs=NA,reps=100,maxlines=6,perc=NA,xlims=NA){

  nm <-deparse(substitute(DLM_data))
  PosMPs<-Can(DLM_data)
  PosMPs<-PosMPs[PosMPs%in%avail("DLM_output")]
  DLM_data@PosMPs<-PosMPs
  if(!is.na(MPs[1]))DLM_data@MPs<-MPs[MPs%in%PosMPs]
  if(is.na(MPs[1]))DLM_data@MPs<-PosMPs
  funcs<-DLM_data@MPs

  if(length(funcs)==0){
    stop("None of the methods 'MPs' are possible given the data available")
  }else{
    TACa<-getTAC(DLM_data,MPs=funcs,reps)
    DLM_data@TAC<-TACa
    return(DLM_data)
    #assign(nm,DLM,envir=.GlobalEnv)
  }

}

getTAC<-function(DLM_data,MPs=NA,reps=100){

  nsims<-length(DLM_data@Mort)
  nMPs<-length(MPs)
  TACa<-array(NA,dim=c(nMPs,reps,nsims))

  if(!sfIsRunning()|(nMPs<8&nsims<8)){
    for(ff in 1:nMPs){
      TACa[ff,,]<-sapply(1:nsims,MPs[ff],DLM_data=DLM_data,reps=reps)
    }
  }else{
    sfExport(list=c("DLM_data"))
    if(nsims<8){
      sfExport(list=c("MPs","reps"))
      for(ss in 1:length(nsims)){
        TACa[,,ss]<-t(sfSapply(1:length(MPs),parallelMPs,DLM_data=DLM_data,reps=reps,MPs=MPs,ss=ss))
      }

    }else{

      for(ff in 1:nMPs){
        TACa[ff,,]<-sfSapply(1:nsims,MPs[ff],DLM_data=DLM_data,reps=reps)
      }
    }
  }
  for(ff in 1:nMPs){
    if(sum(is.na(TACa[ff,,]))>sum(!is.na(TACa[ff,,]))){  # only plot if there are sufficient non-NA TAC samples
      print(paste("Method ",MPs[ff]," produced greater than 50% NA values",sep=""))
    }
  }
  TACa
}

parallelMPs<-function(x,DLM_data,reps,MPs,ss) sapply(ss,MPs[x],DLM_data,reps=reps)

# Create a plot method for DLM data objects
setMethod("plot",
  signature(x = "DLM_data"),
  function(x,funcs=NA,maxlines=6,perc=0.5,xlims=NA){
    
    DLM_data<-x
    cols<-rep(c('black','red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey'),4)
    ltys<-rep(1:4,each=13)
    
    
    if(is.na(funcs[1]))funcs<-DLM_data@MPs

    nMPs<-length(funcs)
    nplots<-ceiling(nMPs/maxlines)
    maxl<-ceiling(nMPs/nplots)
    mbyp <- split(1:nMPs, ceiling(1:nMPs/maxl))   # assign methods to plots

    if(is.na(xlims[1])|length(xlims)!=2){
      xlims<-quantile(DLM_data@TAC,c(0.005,0.95),na.rm=T)
      if(xlims[1]<0)xlims[1]<-0
    }
    if(!NAor0(DLM_data@Ref)){
      if(xlims[1]>DLM_data@Ref)xlims[1]<-max(0,0.98*DLM_data@Ref)
      if(xlims[2]<DLM_data@Ref)xlims[2]<-1.02*DLM_data@Ref
    }
    ylims<-c(0,1)

    #for(m in 1:nMPs){
     # if(sum(!is.na(DLM_data@TAC[m,,1]))>2){
       # dens<-density(DLM_data@TAC[m,,1],na.rm=T)
        #print(quantile(dens$y,0.99,na.rm=T))
      #  if(quantile(dens$y,0.9,na.rm=T)>ylims[2])ylims[2]<-quantile(dens$y,0.90,na.rm=T)
      #}
    #}

    #dev.new2(width=10,height=0.5+7*nplots)
    par(mfrow=c(ceiling(nplots/2),2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.35,0.35,0.35,0.05))

    for(p in 1:nplots){
      m<-mbyp[[p]][1]
      plot(NA,NA,xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col="white",lwd=3,type="l")
      abline(h=0)
      if(!NAor0(DLM_data@Ref)){
        abline(v=DLM_data@Ref,col="light grey",lwd=2)
        if(!NAor0(DLM_data@Ref_type[1]))legend('right',DLM_data@Ref_type,text.col="grey",bty='n')
      }
      #plot(density(DLM@TAC[m,,1],from=0,na.rm=T),xlim=xlims,ylim=ylims,main="",xlab="",ylab="",col=coly[m],lty=ltyy[m],type="l")

      if(!is.na(perc[1]))abline(v=quantile(DLM_data@TAC[m,,1],p=perc,na.rm=T),col=cols[m],lty=ltys[m])
      #if(length(mbyp[[p]])>0){
        for(ll in 1:length(mbyp[[p]])){
          m<-mbyp[[p]][ll]
          if(sum(!is.na(DLM_data@TAC[m,,1]))>10){  # only plot if there are sufficient non-NA TAC samples
            x<-density(DLM_data@TAC[m,,1],from=0,na.rm=T)$x
            y<-density(DLM_data@TAC[m,,1],from=0,na.rm=T)$y
            y<-y/max(y)
            lines(x,y,col=cols[ll])
          }else{
            print(paste("Method ",funcs[m]," produced too many NA TAC values for plotting densities",sep=""))
          }
          if(!is.na(perc[1]))abline(v=quantile(DLM_data@TAC[m,,1],p=perc,na.rm=T),col=cols[ll],lty=2)
        }
      #}
      cind<-1:length(mbyp[[p]])
      legend('topright',funcs[mbyp[[p]]],text.col=cols[cind],col=cols[cind],lty=1,bty='n',cex=0.75)
    }

    mtext(paste("TAC (",DLM_data@Units,")",sep=""),1,outer=T,line=0.5)
    mtext(paste("Standardized relative frequency",sep=""),2,outer=T,line=0.5)
    mtext(paste("TAC calculation for ",DLM_data@Name,sep=""),3,outer=T,line=0.5)
})

condmet<-function(vec)TRUE%in%vec

Sense<-function(DLM_data,MP,nsense=6,reps=100,perc=c(0.05,0.5,0.95),ploty=T){

  DLM_data2<-DLM_data
  nm <-deparse(substitute(DLM_data2))
  refTAC<-quantile(getTAC(DLM_data2,MP,reps),perc,na.rm=T)
 
  DLM_data<-DLM_data2
  reqs<-Required(MP)#read.csv(paste(getwd(),"/Data/Data requirements.csv",sep=""),header=T)
  ind<-(1:nrow(reqs))[reqs[,match(MP,names(reqs))]=="Y"]
  #for(i in 1:length(reqs))
  
 
  slotsCV<-slotNames('DLM_data')[grep("CV_",slotNames('DLM_data'))]  
  slots<-rep("",length(slotsCV))
  for(i in 1:length(slotsCV))slots[i]<-substr(slotsCV[i],4,nchar(slotsCV[i]))
   
  ind<-slots%in%unlist(strsplit(reqs[2],", "))
  slots<-slots[ind]
  slotsCV<-slotsCV[ind]
  sname<-slots
  nslots<-length(slots)

  nrep<-nslots*nsense
  DLM_data<-replic8(DLM_data,nrep)
  pss<-seq(0,1,length.out=nsense+2)[2:(nsense+1)]
  vals<-array(NA,dim=c(nslots,nsense))

  for(i in 1:nslots){
    ind<-(((i-1)*nsense+1):(i*nsense))
    mn<-attr(DLM_data,slots[i])[1]
    cv<-attr(DLM_data,slotsCV[i])[1]*2 # twice the CV of the variable specified in the DLM object
    if(class(attr(DLM_data,slots[i]))=='numeric'){
      if(mn>0){
        attr(DLM_data,slots[i])[ind]<-qlnorm(pss,mconv(mn,cv*mn),sdconv(mn,cv*mn))
        vals[i,]<-qlnorm(pss,mconv(mn,cv*mn),sdconv(mn,cv*mn))
      }else{
        attr(DLM_data,slots[i])[ind]<--qlnorm(pss,mconv(-mn,cv*-mn),sdconv(-mn,cv*-mn))
        vals[i,]<--qlnorm(pss,mconv(-mn,cv*-mn),sdconv(-mn,cv*-mn))
      }
    }else{
      cv<-attr(DLM_data,slotsCV[i])[1]
      attr(DLM_data,slots[i])[ind,]<-attr(DLM_data,slots[i])[ind,]*qlnorm(pss,mconv(1,cv),sdconv(1,cv))
      vals[i,]<-qlnorm(pss,mconv(1,cv),sdconv(1,cv))
    }
  }

  TACa<-getTAC(DLM_data,MPs=MP,reps=reps)
  TACa<-apply(TACa,3,quantile,p=perc,na.rm=T)
  LB<-((1:nslots)-1)*4+1
  UB<-(1:nslots)*4
  sense<-matrix(data=NA,nrow=4*nslots,ncol=nsense+1)
  
  for(i in 1:nslots){
    ind<-((i-1)*nsense+1):(i*nsense)
    dat<-TACa[,ind]

    sense[LB[i],2:(nsense+1)]<-vals[i,]
    sense[(LB[i]+1):UB[i],2:(nsense+1)]<-dat
    sense[LB[i],1]<-slots[i]
    sense[(LB[i]+1):UB[i],1]<-perc
  }
  
  DLM_data2@Sense<-sense

  if(ploty){
   ylimy<-range(TACa)
   #dev.new2(width=10,height=0.5+3*ceiling(nslots/2))
   par(mfrow=c(ceiling(nslots/2),2),mai=c(0.4,0.4,0.01,0.01),omi=c(0.4,0.4,0.4,0.01))
   for(i in 1:nslots){
     ind<-(((i-1)*nsense+1):(i*nsense))
     dat<-TACa[,ind]
     xlimy<-range(vals[i,])
     plot(xlimy,rep(refTAC[2],2),ylim=ylimy,xlim=xlimy,type='l',col="#99999960",main="",xlab="",ylab="")
     abline(h=refTAC[c(1,3)],col="#99999960",lty=2)
     abline(v=slot(DLM_data2,slots[i]),col="#99999960",lty=2)
     lines(vals[i,],dat[2,],col="red",lwd=1.5)
     lines(vals[i,],dat[1,],col="red",lty=2,lwd=1.5)
     lines(vals[i,],dat[3,],col="red",lty=2,lwd=1.5)
     legend('top',legend=sname[i],text.col='blue',bty='n')
   }

   mtext(paste("Output control (",DLM_data@Units,")",sep=""),2,outer=T,line=0.5)
   mtext("Parameter / variable input level",1,outer=T,line=0.5)
   mtext(paste("Sensitivity analysis for ",DLM_data@Name,": ",MP,sep=""),3,outer=T,line=0.5)
  }
  #assign(nm,DLM2,envir=.GlobalEnv)
  DLM_data2
}

replic8<-function(DLM_data,nrep){

  slotnam<-slotNames(DLM_data)
  slotnam<-slotnam[slotnam!="Ref"&slotnam!="OM"&slotnam!="MaxAge"&slotnam!="CAL_bins"&slotnam!="Year"]
  
  for(sl in 1:length(slotnam)){
    slt<-attr(DLM_data,slotnam[sl])
    if(class(slt)=='matrix'){
      attr(DLM_data,slotnam[sl])<-matrix(rep(slt,each=nrep),nrow=nrep,ncol=ncol(slt))
    }else if(class(slt)=='numeric'){
      attr(DLM_data,slotnam[sl])<-rep(slt,nrep)
    }else if(class(slt)=='array'){
      attr(DLM_data,slotnam[sl])<-array(rep(slt,each=nrep),dim=c(nrep,dim(slt)[2:3]))
    }
  }
  DLM_data
}


setMethod("summary",
          signature(object = "DLM_data"),
          function(object){
  
  scols<-c('red','green','blue','orange','brown','purple','dark grey','violet','dark red','pink','dark blue','grey')
  
  #dev.new2(width=8,height=4.5)
  par(mai=c(0.35,0.9,0.2,0.01),c(0.3,0,0,0))
  layout(matrix(c(1,2,1,2,1,2,3,3,3,3),nrow=2))
  plot(object@Year,object@Cat[1,],col="blue",type="l",xlab="Year",ylab=paste("Catch (",object@Units,")",sep=""),ylim=c(0,max(object@Cat[1,],na.rm=T)))
  plot(object@Year,object@Ind[1,],col="orange",type="l",xlab="Year",ylab="Relative abundance",ylim=c(0,max(object@Ind[1,],na.rm=T)))
  
  slots<-c("Dep","Mort","FMSY_M","Dt","BMSY_B0","vbK")
  namey<-c("Stock depletion", "Natural Mortality rate","Ratio of FMSY to M","Depletion over time t","BMSY relative to unfished","Von B. k parameter")
  slotsCV<-c("CV_Dep","CV_Mort","CV_FMSY_M","CV_Dt","CV_BMSY_B0","CV_vbK")
  
  ind<-rep(TRUE,length(slotsCV))
  for(i in 1:length(slotsCV))if(NAor0(attr(object,slots[i]))|NAor0(attr(object,slotsCV[i])))ind[i]<-FALSE
  slots<-slots[ind]
  slotsCV<-slotsCV[ind]
  nrep<-150
  xstore<-array(NA,c(length(slots),nrep))
  ystore<-array(NA,c(length(slots),nrep))
  
 
  for(i in 1:length(slots)){
    mu<-attr(object,slots[i])
    cv<-attr(object,slotsCV[i])
    xstore[i,]<-qlnorm(seq(0,1,length.out=nrep),mconv(mu,cv),sdconv(mu,cv))
    ystore[i,]<-dlnorm(xstore[i,],mconv(mu,cv),sdconv(mu,cv))
  }
  
  plot(xstore[1,],ystore[1,],type="l",xlim=c(0,1.2),ylim=c(0,quantile(ystore,0.97)),xlab="",ylab="Relative frequency",col=scols[1])
  if(length(slots)>1){
    for(i in 2:length(slots)) lines(xstore[i,],ystore[i,],col=scols[i])
  }
  legend('topright',legend=namey[ind],text.col=scols[1:length(slots)],bty='n')
  mtext(paste("Data summary for",deparse(substitute(DLM_data)),sep=" "),3,font=2,line=0.25,outer=T)

})

DLMDataDir<-function(stock=NA){
  if(is.na(stock)){
    paste(searchpaths()[match("package:DLMtool",search())],"/",sep="")
  }else{
    paste(searchpaths()[match("package:DLMtool",search())],"/",stock,".csv",sep="")
  }
}

OneRep<-function(DLM_data){
  DLM_data@CV_Cat=DLM_data@CV_Dt=DLM_data@CV_AvC=DLM_data@CV_Ind=DLM_data@CV_Mort=DLM_data@CV_FMSY_M=DLM_data@CV_BMSY_B0=DLM_data@CV_Cref=DLM_data@CV_Bref=DLM_data@CV_Iref=DLM_data@CV_Rec=DLM_data@CV_Dep=DLM_data@CV_Abun=DLM_data@CV_vbK=DLM_data@CV_vbLinf=DLM_data@CV_vbt0=DLM_data@CV_AM=DLM_data@CV_LFC=DLM_data@CV_LFS=DLM_data@CV_wla=DLM_data@CV_wlb=DLM_data@CV_steep=DLM_data@sigmaL=tiny
  DLM_data
}


#dev.new2 <- function(width = 7, height = 7){
#  platform <- sessionInfo()$platform
#  if (grepl("linux",platform)) {
#    x11(width=width, height=height) 
#  }else if (grepl("pc",platform)) {
#    windows(width=width, height=height)
#  }else if (grepl("apple", platform)) {
#    quartz(width=width, height=height)
#  }
#}

#source('Source code/Fordesktop.R')
