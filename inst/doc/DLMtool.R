## ----Prerequisites_library,echo=TRUE,eval=TRUE---------------------------
library(DLMtool)

## ----Prerequisites_sendall,echo=TRUE-------------------------------------
for(i in 1:length(DLMdat))assign(DLMdat[[i]]@Name,DLMdat[[i]])

## ----Prerequisites_sfinit,echo=TRUE,warning=F----------------------------
sfInit(parallel=T,cpus=2) 

## ----Prerequisites_exportall,echo=TRUE-----------------------------------
sfExportAll()            

## ----Prerequisites_setseed,echo=TRUE-------------------------------------
set.seed(1)            

## ----Demo_operating_model_1,echo=TRUE------------------------------------
OM<-new('OM',                                        
        Rockfish,                                    
        Generic_fleet,                              
        Imprecise_Biased)                            

## ----Demo_operating_model_2,echo=TRUE------------------------------------
slotNames(OM)
class?OM

## ----Demo_methods_1,echo=TRUE--------------------------------------------
avail('DLM quota')
methods<-c("Fratio",                                  
           "DCAC",                          
           "Fdem",    
           "DD")                                      

## ----Demo_methods_2,echo=TRUE--------------------------------------------
?Fratio
?DBSRA

## ----Demo_methods_3,echo=TRUE--------------------------------------------
Fratio

## ----Demo_MSE_1,echo=TRUE------------------------------------------------
RockMSE<-runMSE(OM,methods,nsim=20,reps=1,proyears=30,interval=5)

## ----Demo_MSE_Tplot_1,echo=TRUE,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
Tplot(RockMSE)

## ----Demo_MSE_plot_2,echo=TRUE,fig.width=9,fig.height=4,out.width='0.9\\linewidth',out.height='0.45\\linewidth'----
plot(RockMSE)

## ----Demo_real_data_slots,echo=TRUE--------------------------------------
avail('DLM')
slotNames(Canary_Rockfish)
class?DLM

## ----Demo_real_data_Can_Cant_Needed,echo=TRUE----------------------------
Can(Canary_Rockfish)
Cant(Canary_Rockfish)
Needed(Canary_Rockfish)

## ----Demo_real_getQuota,echo=TRUE----------------------------------------
RockReal<-getQuota(Canary_Rockfish)

## ----Demo_real_pq,echo=TRUE,fig.height=10,fig.width=3,out.width='0.3\\linewidth',out.height='1\\linewidth'----
plot(RockReal)

## ----Demo_real_sense,echo=TRUE,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
RockReal<-Sense(RockReal,"DCAC")

## ----Full_MSE_define_stock,echo=TRUE-------------------------------------
avail('Stock')
ourstock<-Snapper

## ----Full_MSE_setup_bio,echo=TRUE----------------------------------------
ourstock@D<-c(0.05,0.3)
ourstock@Frac_area_1<-c(0.05,0.15)
ourstock@Prob_staying<-c(0.4,0.99)

## ----Full_MSE_setup_fleet,echo=T-----------------------------------------
ourfleet<-Generic_FlatE
ourfleet@Vmaxage<-c(0.5,1)
ourfleet@Spat_targ<-c(1,1.5)

## ----Full_MSE_create_OM,echo=T-------------------------------------------
ourOM<-new('OM',ourstock,ourfleet,Imprecise_Biased)

## ----Full_MSE_runMSE,our_MSE_1,echo=TRUE---------------------------------
ourMSE<-runMSE(ourOM,proyears=20,interval=5,nsim=20,reps=1)

## ----Full_MSE_Tplot,echo=TRUE,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
Tplot(ourMSE)

## ----Full_MSE_subset_targ,echo=TRUE--------------------------------------
Results<-summary(ourMSE) 
Results
Targetted<-subset(Results, Results$Yield>20 & Results$POF<50 & Results$P10<10)
Targetted

## ----Full_MSE_runMSE_2,echo=T--------------------------------------------
ourMSE2<-runMSE(ourOM,Targetted$Method,proyears=20,interval=5,nsim=40,reps=1)

## ----Full_MSE_Pplot,echo=T,fig.width=12,fig.height=11,out.width='0.9\\linewidth',out.height='0.85\\linewidth'----
Pplot(ourMSE2)

## ----Full_MSE_Kplot,echo=T,out.width='0.6\\linewidth',out.height='0.6\\linewidth'----
Kplot(ourMSE2)

## ----Full_MSE_Tplot2,echo=T,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
Tplot(ourMSE2)

## ----Full_MSE_realdata_summary,echo=T,fig.width=7,fig.height=3.5,out.width='0.8\\linewidth',out.height='0.4\\linewidth'----
summary(ourReefFish)

## ----Full_MSE_realdata_getQuota,echo=T-----------------------------------
ourReefFish<-getQuota(ourReefFish)

## ----Full_MSE_realdata_plot_quota,echo=T,fig.width=3,fig.height=8,out.width='0.3\\linewidth',out.height='0.80\\linewidth'----
plot(ourReefFish)

## ----Full_MSE_realdata_Bt_bias,echo=T------------------------------------
ourOM@Btcv
ourOM@Btbiascv

## ----Full_MSE_realdata_Tplot_2_refresh,echo=T,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
Tplot(ourMSE2)

## ----Full_MSE_realdata_sense_DD,echo=T,out.width='0.5\\linewidth',out.height='0.7\\linewidth'----
ourReefFish<-Sense(ourReefFish,'DD')

## ----New_methods_basic_run,echo=T----------------------------------------
sapply(1,Fdem_CC,Red_snapper,reps=5)

## ----New_methods_AvC,echo=T----------------------------------------------
AvC<-function(x,DLM,reps)rlnorm(reps,log(mean(DLM@Cat[x,],na.rm=T)),0.1) 

## ----New_methods_AvC_export,echo=T---------------------------------------
class(AvC)<-"DLM quota"
environment(AvC) <- asNamespace('DLMtool')
sfExport("AvC")

## ----New_methods_THC,echo=T----------------------------------------------
THC<-function(x,DLM,reps){
  rlnorm(reps,log(DLM@Cat[x,order(DLM@Cat[x,],decreasing=T)[3]]),0.1)
}
class(THC)<-"DLM quota"
environment(THC) <- asNamespace('DLMtool')
sfExport("THC")

## ----New_methods_agelim5,echo=T------------------------------------------
agelim5<-function(x,DLM)c(rep(0,4),rep(1,DLM@MaxAge-4)) 
class(agelim5)<-"DLM size"
environment(agelim5) <- asNamespace('DLMtool')
sfExport("agelim5")

## ----New_methods_area1_50,echo=T-----------------------------------------
area1_50<-function(x,DLM)c(0.5,1) 
class(area1_50)<-"DLM space"
environment(area1_50) <- asNamespace('DLMtool')
sfExport("area1_50")

## ----New_methods_MSE,echo=T----------------------------------------------
new_methods<-c("AvC","THC","agelim5","area1_50")
OM<-new('OM',Porgy, Generic_IncE, Imprecise_Unbiased)
PorgMSE<-runMSE(OM,new_methods,maxF=1,nsim=20,reps=1,proyears=20,interval=5)     

## ----New_methods_MSE_Tplot,echo=T,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
Tplot(PorgMSE)                                        

## ----New_methods_MSE_alt_dep,echo=T--------------------------------------
OM@D
OM@D<-c(0.05,0.3)
PorgMSE2<-runMSE(OM,new_methods,maxF=1,nsim=20,reps=1,proyears=20,interval=5)     

## ----New_methods_MSE_alt_dep_Tplot,echo=T,out.width='0.5\\linewidth',out.height='0.5\\linewidth'----
Tplot(PorgMSE2)

## ----Real_data_slotNames,echo=T------------------------------------------
slotNames('DLM')

## ----Real_data_DLMDataDir,echo=T-----------------------------------------
DLMDataDir()

## ----Real_data_Madeup,echo=T---------------------------------------------
Madeup<-new('DLM')                                  #  Create a blank DLM object
Madeup@Name<-'Test'                                 #  Name it
Madeup@Cat<-matrix(20:11*rlnorm(10,0,0.2),nrow=1)   #  Generate fake catch data
Madeup@Units<-"Million metric tonnes"               #  State units of catch
Madeup@AvC<-mean(Madeup@Cat)                        #  Average catches for time t (DCAC)
Madeup@t<-ncol(Madeup@Cat)                          #  No. yrs for Av. catch (DCAC)
Madeup@Dt<-0.5                                      #  Depletion over time t (DCAC)
Madeup@Dep<-0.5                                     #  Depletion relative to unfished 
Madeup@Mort<-0.1                                    #  Natural mortality rate
Madeup@Abun<-200                                    #  Current abundance
Madeup@FMSY_M<-0.75                                 #  Ratio of FMSY/M
Madeup@AM<-3.5                                      #  Age at maturity
Madeup@BMSY_B0<-0.35                                #  BMSY relative to unfished

## ----Real_data_summary,echo=T,fig.width=7,fig.height=3.5,out.width='0.8\\linewidth',out.height='0.4\\linewidth'----
summary(Atlantic_mackerel)

## ----Real_data_CCNR,echo=T-----------------------------------------------
Can(Atlantic_mackerel)
Cant(Atlantic_mackerel)
Needed(Atlantic_mackerel)

## ----Real_data_getQuota,echo=T-------------------------------------------
Atlantic_mackerel<-getQuota(Atlantic_mackerel,reps=48)

## ----Real_data_plot_Quota,echo=T,fig.width=3,fig.height=8,out.width='0.3\\linewidth',out.height='0.80\\linewidth'----
plot(Atlantic_mackerel)

## ----MPA_spec,echo=T-----------------------------------------------------
Rock<-Rockfish
Rock@Prob_staying<-c(0.9,0.999)
Rock@Frac_area_1<-c(0.05,0.5)

## ----MPA_run,echo=T------------------------------------------------------
RockMPA<-runMSE(new('OM',Rock,Generic_fleet,Generic_obs),"area1MPA",nsim=20)
summary(RockMPA)

## ----MPA_plot,echo=TRUE,fig.width=9,fig.height=5,out.width='0.9\\linewidth',out.height='0.45\\linewidth'----
plot(RockMPA)

## ----sfstop,echo=TRUE,include=FALSE,cache=FALSE--------------------------
sfStop()                            

