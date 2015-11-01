## ----Prereq_library,echo=TRUE,eval=TRUE----------------------------------
library(DLMtool)

## ----Prereq_assigndat,echo=TRUE,warning=F--------------------------------
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]])

## ----Prerequisites_sfinit,echo=TRUE,warning=F----------------------------
#sfInit(parallel=T,cpus=2) 

## ----Prerequisites_sfExport,echo=TRUE,warning=F--------------------------
#sfExportAll()

## ----Prereq_random_seed,echo=TRUE,warning=F------------------------------
set.seed(1) 

## ----Demo_operating_model_1,echo=TRUE------------------------------------
OM<-new('OM',                                        
        Blue_shark,                                    
        Generic_fleet,                              
        Imprecise_Biased)                            

## ----Demo_operating_model_2,echo=TRUE------------------------------------
slotNames(OM)
class?OM

## ----Demo_methods_1,echo=TRUE--------------------------------------------
avail('DLM_output')
MPs<-c("Fratio",                                  
           "DCAC",                          
           "Fdem",    
           "DD")                                      

## ----Demo_MPs_2,echo=TRUE------------------------------------------------
?Fratio
?DBSRA

## ----Demo_MPs_3,echo=TRUE------------------------------------------------
Fratio

## ----Demo_MSE_1,echo=TRUE------------------------------------------------
SnapMSE<-runMSE(OM,MPs,nsim=16,reps=1,proyears=30,interval=5)

## ----Demo_MSE_plot_1,echo=TRUE,fig.width=9,fig.height=6,out.width='0.9\\linewidth',out.height='0.6\\linewidth'----
plot(SnapMSE)

## ----Demo_real_data_slots,echo=TRUE--------------------------------------
avail('DLM_data')
slotNames(China_rockfish)
class?DLM_data

## ----Demo_real_data_Can_Cant_Needed,echo=TRUE----------------------------
Can(China_rockfish)
Cant(China_rockfish)
Needed(China_rockfish)

## ----Demo_real_getTAC,echo=TRUE------------------------------------------
RockReal<-TAC(China_rockfish)

## ----Demo_real_pq,echo=TRUE,out.width='0.7\\linewidth',out.height='0.5\\linewidth'----
plot(RockReal)

## ----Demo_real_sense,echo=TRUE,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
RockReal<-Sense(RockReal,"DCAC")

## ----Full_MSE_define_stock,echo=TRUE-------------------------------------
avail('Stock')
ourstock<-Snapper

## ----Full_MSE_setup_bio,echo=TRUE----------------------------------------
ourstock@M<-c(0.2,0.25)
ourstock@maxage<-18
ourstock@D<-c(0.05,0.3)
ourstock@Frac_area_1<-c(0.05,0.15)
ourstock@Prob_staying<-c(0.4,0.99)

## ----Full_MSE_setup_fleet,echo=T-----------------------------------------
ourfleet<-Generic_FlatE
ourfleet@Vmaxlen<-c(0.5,1)

## ----Full_MSE_create_OM,echo=T-------------------------------------------
ourOM<-new('OM',ourstock,ourfleet,Imprecise_Biased)

## ----Full_MSE_runMSE,our_MSE_1,echo=TRUE---------------------------------
ourMSE<-runMSE(ourOM,proyears=20,interval=5,nsim=16,reps=1)

## ----Full_MSE_Tplot,echo=TRUE,out.width='0.8\\linewidth',out.height='0.8\\linewidth'----
Tplot(ourMSE)

## ----Full_MSE_subset_targ,echo=TRUE--------------------------------------
Results<-summary(ourMSE) 
Results
Targetted<-subset(Results, Results$Yield>50 & Results$POF<50 & Results$P50<20)
Targetted

## ----Full_MSE_runMSE_2,echo=T--------------------------------------------
TargMP<-Targetted$MP[grep("FMSYref",Targetted$MP,invert=T)]
ourMSE2<-runMSE(ourOM,TargMP,proyears=20,interval=5,nsim=32,reps=1)

## ----Full_MSE_CheckConverge,echo=T,fig.width=12,fig.height=8,out.width='0.9\\linewidth',out.height='0.65\\linewidth'----
CheckConverg(ourMSE2)

## ----Full_MSE_Pplot,echo=T,fig.width=12,fig.height=8,out.width='0.9\\linewidth',out.height='0.65\\linewidth'----
Pplot(ourMSE2)

## ----Full_MSE_Kplot,echo=T,out.width='0.8\\linewidth',out.height='0.6\\linewidth'----
Kplot(ourMSE2)

## ----Full_MSE_Tplot2,echo=T,out.width='0.7\\linewidth',out.height='0.35\\linewidth'----
Tplot2(ourMSE2)

## ----Full_MSE_VOI,echo=T,fig.width=12,fig.height=11,out.width='0.9\\linewidth',out.height='0.85\\linewidth'----
VOI(ourMSE2)

## ----Full_MSE_realdata_summary,echo=T,fig.width=7,fig.height=3.5,out.width='0.8\\linewidth',out.height='0.4\\linewidth'----
summary(ourReefFish)

## ----Full_MSE_realdata_getTAC,echo=T-------------------------------------
ourReefFish<-TAC(ourReefFish)

## ----Full_MSE_realdata_plot_TAC,echo=T,fig.width=6,fig.height=8,out.width='0.7\\linewidth',out.height='0.80\\linewidth'----
plot(ourReefFish)

## ----Full_MSE_realdata_Bt_bias,echo=T------------------------------------
ourOM@Btcv
ourOM@Btbias

## ----Full_MSE_realdata_Tplot_2_refresh,echo=T,out.width='0.7\\linewidth',out.height='0.35\\linewidth'----
Tplot2(ourMSE2)

## ----Full_MSE_realdata_sense_DD,echo=T,out.width='0.5\\linewidth',out.height='0.7\\linewidth'----
ourReefFish<-Sense(ourReefFish,'DD')

## ----New_methods_basic_run,echo=T----------------------------------------
sapply(1,Fdem_CC,Red_snapper,reps=5)

## ----New_MPs_AvC,echo=T--------------------------------------------------
AvC<-function(x,DLM_data,reps)rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),0.1) 

## ----New_MPs_AvC_export,echo=T-------------------------------------------
class(AvC)<-"DLM_output"
environment(AvC) <- asNamespace('DLMtool')
#sfExport("AvC")

## ----New_MPs_THC,echo=T--------------------------------------------------
THC<-function(x,DLM_data,reps){
  rlnorm(reps,log(DLM_data@Cat[x,order(DLM_data@Cat[x,],decreasing=T)[3]]),0.1)
}
class(THC)<-"DLM_output"
environment(THC) <- asNamespace('DLMtool')
#sfExport("THC")

## ----New_MPs_agelim5,echo=T----------------------------------------------
agelim5<-function(x,DLM_data){
  Allocate<-1 # Fraction of effort reallocated to open area
  Effort<-1  # Fraction of effort in last historical year
  Spatial<-c(1,1) # Fraction of effort found in each area 
  Vuln<-c(rep(0,4),rep(1,DLM_data@MaxAge-4)) # Age vulnerability 
  c(Allocate, Effort, Spatial, Vuln) # Input controls stitched togther
}
class(agelim5)<-"DLM_input"
environment(agelim5) <- asNamespace('DLMtool')
#sfExport("agelim5")

## ----New_MPs_area1_50,echo=T---------------------------------------------
area1_50<-function(x,DLM_data){ 
  Allocate<-0 # Fraction of effort reallocated to open area
  Effort<-1  # Fraction of effort in last historical year
  Spatial<-c(0.5,1) # Fraction of effort found in each area
  Vuln<-rep(NA,DLM_data@MaxAge) # Age vulnerability is not specified   
  c(Allocate, Effort, Spatial, Vuln) # Input controls stitched togther
}
class(area1_50)<-"DLM_input"
environment(area1_50) <- asNamespace('DLMtool')
#sfExport("area1_50")

## ----New_MPs_MSE,echo=T--------------------------------------------------
new_MPs<-c("AvC","THC","agelim5","area1_50")
OM<-new('OM',Porgy, Generic_IncE, Imprecise_Unbiased)
PorgMSE<-runMSE(OM,new_MPs,maxF=1,nsim=20,reps=1,proyears=20,interval=5)     

## ----New_MPs_MSE_Tplot,echo=T,out.width='0.7\\linewidth',out.height='0.7\\linewidth'----
Tplot(PorgMSE)                                        

## ----New_MPs_MSE_alt_dep,echo=T------------------------------------------
OM@D
OM@D<-c(0.05,0.3)
PorgMSE2<-runMSE(OM,new_MPs,maxF=1,nsim=16,reps=1,proyears=20,interval=5)     

## ----New_MPs_MSE_alt_dep_Tplot,echo=T,out.width='0.7\\linewidth',out.height='0.7\\linewidth'----
Tplot(PorgMSE2)

## ----Real_data_slotNames,echo=T------------------------------------------
slotNames('DLM_data')

## ----Real_data_DLMDataDir,echo=T-----------------------------------------
DLMDataDir()

## ----Real_data_Madeup,echo=T---------------------------------------------
Madeup<-new('DLM_data')                             #  Create a blank DLM object
Madeup@Name<-'Test'                                 #  Name it
Madeup@Cat<-matrix(20:11*rlnorm(10,0,0.2),nrow=1)   #  Generate fake catch data
Madeup@Units<-"Million metric tonnes"               #  State units of catch
Madeup@AvC<-mean(Madeup@Cat)                        #  Average catches for time t (DCAC)
Madeup@t<-ncol(Madeup@Cat)                          #  No. yrs for Av. catch (DCAC)
Madeup@Dt<-0.5                                      #  Depletion over time t (DCAC)
Madeup@Dep<-0.5                                     #  Depletion relative to unfished 
Madeup@vbK<-0.2                                     #  VB maximum growth rate
Madeup@vbt0<-(-0.5)                                 #  VB theoretical age at zero length
Madeup@vbLinf<-200                                  #  VB maximum length
Madeup@Mort<-0.1                                    #  Natural mortality rate
Madeup@Abun<-200                                    #  Current abundance
Madeup@FMSY_M<-0.75                                 #  Ratio of FMSY/M
Madeup@L50<-100                                     #  Length at 50% maturity
Madeup@L95<-120                                     #  Length at 95% maturity
Madeup@BMSY_B0<-0.35                                #  BMSY relative to unfished

## ----Real_data_summary,echo=T,fig.width=7,fig.height=3.5,out.width='0.8\\linewidth',out.height='0.4\\linewidth'----
summary(Atlantic_mackerel)

## ----Real_data_CCNR,echo=T-----------------------------------------------
Can(Atlantic_mackerel)
Cant(Atlantic_mackerel)
Needed(Atlantic_mackerel)

## ----Real_data_getTAC,echo=T---------------------------------------------
Atlantic_mackerel<-TAC(Atlantic_mackerel,reps=48)

## ----Real_data_plot_TAC,echo=T,fig.width=6,fig.height=8,out.width='0.6\\linewidth',out.height='0.80\\linewidth'----
plot(Atlantic_mackerel)

## ----sfstop,echo=TRUE,include=FALSE,cache=FALSE--------------------------
#sfStop()                            

