# Input MPs

#' A data-limited method in which fishing vulnerability is set according to the
#' maturity curve
#' 
#' An example of the implementation of input controls in the DLM toolkit, where
#' selectivity-at-length is set equivalent to maturity-at-length
#' 
#' 
#' @usage matlenlim(x, Data, ...)
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return A vector of input control recommendations, with values for length at
#' first capture and full selection
#' @author T. Carruthers
#' @references Made-up for this package
#' @export matlenlim
matlenlim <- function(x, Data, ...) {
  # Length at maturity is knife-edge vulnerability
  dependencies = "Data@LFC, Data@LFS"
  Allocate <- 1
  Effort <- 1
  Spatial <- c(1, 1)
  
  newLFC <- Data@L50[x] * 0.95
  newLFS <- Data@L50[x]
  Vuln <- c(newLFC, newLFS)
  c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim) <- "Input"



#' A data-limited method in which fishing vulnerability is set slightly higher
#' than the maturity curve
#' 
#' An example of the implementation of input controls in the DLM toolkit, where
#' selectivity-at-length is set slightly higher than the maturity-at-length
#' 
#' 
#' @usage matlenlim2(x, Data, ...)
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return A vector of input control recommendations, with values for length at
#' first capture and full selection
#' @author A. Hordyk
#' @references Made-up for this package
#' @export matlenlim2
matlenlim2 <- function(x, Data, ...) {
  # Knife-edge vulnerability slightly higher than length at maturity
  dependencies = "Data@LFC, Data@LFS"
  Allocate <- 1
  Effort <- 1
  Spatial <- c(1, 1)
  newLFS <- 1.1 * Data@L50[x]
  newLFC <- 0.95 * newLFS
  Vuln <- c(newLFC, newLFS)
  c(Allocate, Effort, Spatial, Vuln)
}
class(matlenlim2) <- "Input"



#' An data-limited method which sets a slot limit
#' 
#' An example of the implementation of input controls in the DLM toolkit, where
#' selectivity-at-length is set using a slot limit; that is, a minimum and
#' maximum legal length.  The maximum limit is set here, quite arbitrarily, as
#' the 75th percentile between the new minimum legal length and the estimated
#' asymptotic length.
#' 
#' 
#' @usage slotlim(x, Data, ...)
#' @param x A position in a data-limited methods object
#' @param Data A data-limited methods object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @return A vector of input control recommendations, with values for length at
#' first capture, full selection, and maximum size limit in the 5th, 6th, and
#' 7th elements of the vector
#' @author A. Hordyk
#' @references Made-up for this package
#' @export slotlim
slotlim <- function(x, Data, ...) {
  # Example of slot limit between 0.95 and 1.25 * L50
  dependencies = "Data@LFC, Data@LFS"
  Allocate <- 1
  Effort <- 1
  Spatial <- c(1, 1)
  
  newLFS <- 14 + 1.1 * Data@L50[x]
  newLFC <- 0.95 * newLFS
  UppLim <- as.numeric(quantile(c(newLFS, Data@vbLinf[x]), 0.75))
  Vuln <- c(newLFC, newLFS, UppLim)
  c(Allocate, Effort, Spatial, Vuln)
}
class(slotlim) <- "Input"



#' An marine reserve in area 1 with full reallocation of fishing effort
#' 
#' A spatial control that prevents fishing in area 1 and reallocates this
#' fishing effort to area 2.
#' 
#' 
#' @usage MRreal(x, Data, ...)
#' @param x A position in data / simulation object DLM
#' @param Data A data limited methods data object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @author T. Carruthers
#' @export MRreal
MRreal <- function(x, Data, ...) {
  # A Marine reserve in area 1 with spatial reallocation of effort
  dependencies = "Data@MaxAge"
  Allocate <- 1
  Effort <- 1
  Spatial <- c(0, 1)
  # Vuln<-rep(NA,Data@MaxAge)
  Vuln <- rep(NA, 2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRreal) <- "Input"



#' An marine reserve in area 1 with no spatial reallocation of fishing effort
#' 
#' A spatial control that prevents fishing in area 1 and does not reallocate
#' this fishing effort to area 2.
#' 
#' 
#' @usage MRnoreal(x, Data, ...)
#' @param x A position in data / simulation object DLM
#' @param Data A data limited methods data object
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @author T. Carruthers
#' @export MRnoreal
MRnoreal <- function(x, Data, ...) {
  # A Marine reserve in area 1 with no spatial reallocation of effort
  dependencies = "Data@MaxAge"
  Allocate <- 0
  Effort <- 1
  Spatial <- c(0, 1)
  # Vuln<-rep(NA,Data@MaxAge)
  Vuln <- rep(NA, 2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(MRnoreal) <- "Input"



#' Fishing at current effort levels
#' 
#' Constant fishing effort set at final year of historical simulations subject
#' to changes in catchability determined by OM@qinc and interannual variability
#' in catchability determined by OM@qcv. This MP is intended to represent a
#' 'status quo' management approach.
#' 
#' 
#' @usage curE(x, Data, ...)
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @note Made up for this package.
#' @author T. Carruthers.
#' @export curE
curE <- function(x, Data, ...) {
  # current effort
  dependencies = "Data@MaxAge"
  Allocate <- 1
  Effort <- 1
  Spatial <- c(1, 1)
  # Vuln<-rep(NA,Data@MaxAge)
  Vuln <- rep(NA, 2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE) <- "Input"



#' Fishing at 75 per cent of current effort levels
#' 
#' Constant fishing effort set at 75 per cent of final year of historical
#' simulations subject to changes in catchability determined by OM@qinc and
#' interannual variability in catchability determined by OM@qcv. This MP is
#' intended to represent a 'status quo' management approach.
#' 
#' 
#' @usage curE75(x, Data, ...)
#' @param x A position in a data-limited methods data object.
#' @param Data A data-limited methods data object.
#' @param ... Optional additional arguments that are ignored. Note arguments
#' \code{reps} or \code{...} are required for all input controls
#' @note Made up for this package.
#' @author T. Carruthers.
#' @export curE75
curE75 <- function(x, Data, ...) {
  # 75% current effort
  dependencies = "Data@MaxAge"
  Allocate <- 1
  Effort <- 0.75
  Spatial <- c(1, 1)
  # Vuln<-rep(NA,Data@MaxAge)
  Vuln <- rep(NA, 2)
  c(Allocate, Effort, Spatial, Vuln)
}
class(curE75) <- "Input"





#' Effort control version of DD - Delay - Difference Stock Assessment with UMSY
#' and MSY leading
#' 
#' A simple delay-difference assessment that estimates and recommends FMSY
#' using a time-series of catches and a relative abundance index.
#' 
#' 
#' @usage DDe(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export DDe
DDe <- function(x, Data, reps = 100) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
    Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
    ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
    Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
    C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
    upper = log(exp(params) * 20), hessian = TRUE)
  
  U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  
  Allocate <- 1
  eff <- exp(opt$par[1])/U_hist[Data@LHYear]
  eff[!is.finite(eff)] <- 0.01
  eff[eff > 1e+05] <- 0.01
  Effort <- max(0.01, eff)
  
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  # Out <- list() Out[[1]] <- out Out[[2]] <- MiscList
  return(out)
  
}
class(DDe) <- "Input"



#' Effort control version of DD - Delay - Difference Stock Assessment with UMSY
#' and MSY leading that fishes at 75 per cent of FMSY
#' 
#' A simple delay-difference assessment that estimates and recommends 75 per
#' cent FMSY using a time-series of catches and a relative abundance index.
#' 
#' 
#' @usage DDe75(x, Data, reps = 100)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export DDe75
DDe75 <- function(x, Data, reps = 100) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
    Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
    ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
    Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
    C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
    upper = log(exp(params) * 20), hessian = TRUE)
  
  U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  
  Allocate <- 1
  eff <- exp(0.75 * opt$par[1])/U_hist[Data@LHYear]
  eff[!is.finite(eff)] <- 0.01
  eff[eff > 1e+05] <- 0.01
  Effort <- max(0.01, eff)
  
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  # Out <- list() Out[[1]] <- out Out[[2]] <- MiscList
  
  return(out)
  
}
class(DDe75) <- "Input"



#' Effort searching version of DD - Delay - Difference Stock Assessment with
#' UMSY and MSY leading that fishes at 75 per cent of FMSY
#' 
#' A simple delay-difference assessment that estimates FMSY using a time-series
#' of catches and a relative abundance index. The MP provides a change in
#' effort in the direction of FMSY up to a maximum change of 10 percent.
#' 
#' 
#' @usage DDes(x, Data, reps = 100, LB=0.9, UB=1.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @note This DD model is observation error only and has does not estimate
#' process error (recruitment deviations). Similar to many other assessment
#' models it depends on a whole host of dubious assumptions such as temporally
#' stationary productivity and proportionality between the abundance index and
#' real abundance. Unsurprisingly the extent to which these assumptions are
#' violated tends to be the biggest driver of performance for this method.
#' @author T. Carruthers
#' @references Method based on equations of Carl Walters (bug him with
#' questions and expect colourful responses)
#' @export DDes
DDes <- function(x, Data, reps = 100, LB = 0.9, UB = 1.1) {
  # for(x in 1:nsim){
  dependencies = "Data@vbLinf, Data@CV_vbLinf, Data@vbK, Data@CV_vbK, Data@vbt0, Data@CV_vbt0, Data@Mort, Data@CV_Mort. Data@wla, Data@ wlb"
  Linfc <- trlnorm(reps, Data@vbLinf[x], Data@CV_vbLinf[x])
  Kc <- trlnorm(reps, Data@vbK[x], Data@CV_vbK[x])
  if (Data@vbt0[x] != 0 & Data@CV_vbt0[x] != tiny) {
    t0c <- -trlnorm(reps, -Data@vbt0[x], Data@CV_vbt0[x])
  } else {
    t0c <- rep(Data@vbt0[x], reps)
  }
  t0c[!is.finite(t0c)] <- 0
  Mdb <- trlnorm(reps, Data@Mort[x], Data@CV_Mort[x])  # CV of 0.5 as in MacCall 2009
  a <- Data@wla[x]
  b <- Data@wlb[x]
  
  Winf = Data@wla[x] * Data@vbLinf[x]^Data@wlb[x]
  age <- 1:Data@MaxAge
  la <- Data@vbLinf[x] * (1 - exp(-Data@vbK[x] * ((age - Data@vbt0[x]))))
  wa <- Data@wla[x] * la^Data@wlb[x]
  a50V <- iVB(Data@vbt0[x], Data@vbK[x], Data@vbLinf[x], 
    Data@L50[x])
  a50V <- max(a50V, 1)
  yind <- (1:length(Data@Cat[x, ]))[!is.na(Data@Cat[x, ] + Data@Ind[x, 
    ])]
  C_hist <- Data@Cat[x, yind]
  E_hist <- C_hist/Data@Ind[x, yind]
  E_hist <- E_hist/mean(E_hist)
  ny_DD <- length(C_hist)
  params <- log(c(Data@Mort[x], mean(C_hist, na.rm = T), Data@Mort[x]))
  k_DD <- ceiling(a50V)  # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD > Data@MaxAge/2] <- ceiling(Data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD <- (wa[k_DD + 2] - Winf)/(wa[k_DD + 1] - Winf)
  Alpha_DD <- Winf * (1 - Rho_DD)
  So_DD <- exp(-Data@Mort[x])  # get So survival rate
  wa_DD <- wa[k_DD]
  UMSYprior <- c(1 - exp(-Data@Mort[x] * 0.5), 0.3)
  opt <- optim(params, DD_R, opty = 1, So_DD = So_DD, Alpha_DD = Alpha_DD, 
    Rho_DD = Rho_DD, ny_DD = ny_DD, k_DD = k_DD, wa_DD = wa_DD, E_hist = E_hist, 
    C_hist = C_hist, UMSYprior = UMSYprior, method = "L-BFGS-B", lower = log(exp(params)/20), 
    upper = log(exp(params) * 20), hessian = TRUE)
  
  U_hist <- 1 - exp(-exp(opt$par[3]) * E_hist)
  fac <- exp(opt$par[1])/U_hist[Data@LHYear]  # ratio of UMSY to reference U
  fac <- fac * (U_hist[Data@LHYear]/U_hist[length(U_hist)])  # ratio of last U to reference U
  
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  Allocate <- 1
  Effort <- max(0.01, Data@MPeff[x] * fac)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  # Out <- list() Out[[1]] <- out Out[[2]] <- MiscList
  
  return(out)
  
}
class(DDes) <- "Input"



#' Effort searching MP aiming for 40 per cent stock depletion
#' 
#' A very simple MP that modifies effort to reach 40 percent stock depletion
#' 
#' 
#' @usage DTe40(x, Data, reps = 100, alpha=0.4, LB=0.9, UB=1.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param alpha The target level of depletion
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @author T. Carruthers
#' @export DTe40
DTe40 <- function(x, Data, reps = 100, alpha = 0.4, LB = 0.9, UB = 1.1) {
  
  dependencies = "Data@Dep"
  
  fac <- Data@Dep[x]/alpha
  
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  Allocate <- 1
  Effort <- max(0.01, Data@MPeff[x] * fac)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  # MiscList<-Effort
  
  # Out <- list() Out[[1]] <- out Out[[2]] <- MiscList
  
  return(out)
  
}
class(DTe40) <- "Input"



#' Effort searching MP aiming for 50 per cent stock depletion
#' 
#' A very simple MP that modifies effort to reach 50 percent stock depletion
#' 
#' 
#' @usage DTe50(x, Data, reps = 100, alpha=0.5, LB=0.9, UB=1.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the TAC recommendation
#' @param alpha The target level of depletion
#' @param LB The lowest permitted factor of previous fishing effort
#' @param UB The highest permitted factor of previous fishing effort
#' @author T. Carruthers
#' @export DTe50
DTe50 <- function(x, Data, reps = 100, alpha = 0.5, LB = 0.9, UB = 1.1) {
  
  dependencies = "Data@Dep"
  
  fac <- Data@Dep[x]/alpha
  if (fac < LB) 
    fac <- LB
  if (fac > UB) 
    fac <- UB
  
  Allocate <- 1
  Effort <- max(Data@MPeff[x] * fac, 0.01)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  
  # MiscList<-Effort
  
  # Out <- list() Out[[1]] <- out Out[[2]] <- MiscList
  
  return(out)
  
}
class(DTe50) <- "Input"



#' A management procedure that incrementally adjusts the TAC according to the
#' mean length of recent catches.
#' 
#' A effort-based version of least biologically precautionary of four adaptive
#' length-based MPs proposed by Geromont and Butterworth 2014. Tested by
#' Carruthers et al. 2015
#' 
#' 
#' @usage LstepCE1(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.05,
#' llim = c(0.96, 0.98, 1.05))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of effort samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of the effort update increment.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the effort.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LstepCE1
LstepCE1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.05, 
  llim = c(0.96, 0.98, 1.05)) {
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  # ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  rat <- Lrecent/Lave
  
  step <- stepsz
  
  if (rat < llim[1]) {
    Effort <- Data@MPeff[x] - 2 * (step * Data@MPeff[x])
  } else if (rat < llim[2]) {
    Effort <- Data@MPeff[x] - (step * Data@MPeff[x])
  } else if (rat > llim[3]) {
    Effort <- Data@MPeff[x] + (step * Data@MPeff[x])
  } else {
    Effort <- Data@MPeff[x]
  }
  
  Allocate <- 1
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(LstepCE1) <- "Input"



#' A management procedure that incrementally adjusts the Effort according to
#' the mean length of recent catches.
#' 
#' A effort-based version of one of the four adaptive length-based MPs proposed
#' by Geromont and Butterworth 2014.
#' 
#' 
#' @usage LstepCE2(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.1,
#' llim = c(0.96, 0.98, 1.05))
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param stepsz Parameter controlling the size of the effort update increment.
#' @param llim A vector of length reference points that determine the
#' conditions for increasing, maintaining or reducing the effort.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LstepCE2
LstepCE2 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, stepsz = 0.1, 
  llim = c(0.96, 0.98, 1.05)) {
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  # ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  rat <- Lrecent/Lave
  step <- stepsz
  
  if (rat < llim[1]) {
    Effort <- Data@MPeff[x] - 2 * (step * Data@MPeff[x])
  } else if (rat < llim[2]) {
    Effort <- Data@MPeff[x] - (step * Data@MPeff[x])
  } else if (rat > llim[3]) {
    Effort <- Data@MPeff[x] + (step * Data@MPeff[x])
  } else {
    Effort <- Data@MPeff[x]
  }
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Allocate <- 1
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(LstepCE2) <- "Input"



#' A management procedure that incrementally adjusts the Effort to reach a
#' target mean length in catches.
#' 
#' A effort based version of the least biologically precautionary of four
#' target length MPs proposed by Geromont and Butterworth 2014.
#' 
#' 
#' @usage LtargetE1(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.05)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LtargetE1
LtargetE1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.05) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  L0 <- 0.9 * Lave
  Ltarget <- xL * Lave
  if (Lrecent > L0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Lrecent - L0)/(Ltarget - 
      L0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Lrecent/L0)^2
  }
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.85] <- 0.85
  Step[Step > 1.15] <- 1.15
  
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(LtargetE1) <- "Input"



#' A management procedure that incrementally adjusts the Effort to reach a
#' target mean length in catches.
#' 
#' A effort based version of the most biologically precautionary of four target
#' length MPs proposed by Geromont and Butterworth 2014.
#' 
#' 
#' @usage LtargetE4(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.15)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param xL Parameter controlling the magnitude of the target mean length of
#' catches relative to average length in catches.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export LtargetE4
LtargetE4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, xL = 1.15) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Lrecent <- mean(Data@ML[ind])
  Lave <- mean(Data@ML[ind3])
  L0 <- 0.9 * Lave
  Ltarget <- xL * Lave
  if (Lrecent > L0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Lrecent - L0)/(Ltarget - 
      L0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Lrecent/L0)^2
  }
  
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.8] <- 0.8
  Step[Step > 1.2] <- 1.2
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(LtargetE4) <- "Input"



#' A management procedure that incrementally adjusts the effort to reach a
#' target CPUE / relative abundance index
#' 
#' An effort-based version of the least biologically precautionary of two
#' index/CPUE target MPs proposed by Geromont and Butterworth 2014. Tested by
#' Carruthers et al. 2015
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage ItargetE1(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @export ItargetE1
ItargetE1 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 1.5) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Irecent - I0)/(Itarget - 
      I0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Irecent/I0)^2
  }
  
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.85] <- 0.85
  Step[Step > 1.15] <- 1.15
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(ItargetE1) <- "Input"



#' A management procedure that incrementally adjusts the Effort to reach a
#' target CPUE / relative abundance index
#' 
#' An effort-based version of the most biologically precautionary of two
#' index/CPUE target MPs proposed by Geromont and Butterworth 2014.
#' 
#' Tested by Carruthers et al. 2015.
#' 
#' @usage ItargetE4(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 2.5)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of samples
#' @param yrsmth Years over which to smooth recent estimates of surplus
#' production
#' @param xx Parameter controlling the fraction of mean catch to start using in
#' first year
#' @param Imulti Parameter controlling how much larger target CPUE / index is
#' compared with recent levels.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @references Carruthers et al. 2015. Performance evaluation of simple
#' management procedures. Fish and Fisheries. In press.
#' 
#' Geromont, H.F., Butterworth, D.S. 2014. Generic management procedures for
#' data-poor fisheries; forecasting with few data. ICES J. Mar. Sci.
#' doi:10.1093/icesjms/fst232
#' @export ItargetE4
ItargetE4 <- function(x, Data, reps = 100, yrsmth = 5, xx = 0, Imulti = 2.5) {
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 5 years
  ylast <- (Data@LHYear - Data@Year[1]) + 1  #last historical year
  ind2 <- ((ylast - (yrsmth - 1)):ylast)  # historical 5 pre-projection years
  ind3 <- ((ylast - (yrsmth * 2 - 1)):ylast)  # historical 10 pre-projection years
  
  Irecent <- mean(Data@Ind[x, ind])
  Iave <- mean(Data@Ind[x, ind3])
  Itarget <- Iave * Imulti
  I0 <- 0.8 * Iave
  if (Irecent > I0) {
    Effort <- 0.5 * Data@MPeff[x] * (1 + ((Irecent - I0)/(Itarget - 
      I0)))
  } else {
    Effort <- 0.5 * Data@MPeff[x] * (Irecent/I0)^2
  }
  Step <- (Effort/Data@MPeff[x])  # step change in effort 
  Step[Step < 0.8] <- 0.8
  Step[Step > 1.2] <- 1.2
  
  Allocate <- 1
  Effort <- Step * Data@MPeff[x]
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(ItargetE4) <- "Input"



#' Index Target Effort-Based 10
#' 
#' An index target MP where the Effort is modified according to current index
#' levels (mean index over last 5 years) relative to a target level. Maximum
#' annual changes are 10 per cent.
#' 
#' 
#' @usage ITe10(x, Data, reps = 100, yrsmth = 5, mc = 0.1)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the Effort among years.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export ITe10
ITe10 <- function(x, Data, reps = 100, yrsmth = 5, mc = 0.1) {
  
  dependencies = "Data@Ind, Data@MPeff, DLMdata@CV_Ind, DLMdata@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) 
    deltaI <- 1 - mc
  if (deltaI > (1 + mc)) 
    deltaI <- 1 + mc
  
  Effort <- Data@MPeff[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  if (reps == 1) 
    Effort <- Data@MPeff[x] * deltaI
  Allocate <- 1
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(ITe10) <- "Input"



#' Index Target Effort-Based 5
#' 
#' An index target MP where the Effort is modified according to current index
#' levels (mean index over last 5 years) relative to a target level. Maximum
#' annual changes are 5 per cent.
#' 
#' 
#' @usage ITe5(x, Data, reps = 100, yrsmth = 5, mc = 0.05)
#' @param x A position in a data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of stochastic samples of the quota recommendation
#' @param yrsmth The number of historical years over which to average the index
#' @param mc The maximum fractional change in the effort among years.
#' @return A numeric vector of input controls
#' @author T. Carruthers
#' @export ITe5
ITe5 <- function(x, Data, reps = 100, yrsmth = 5, mc = 0.05) {
  
  dependencies = "Data@Ind, Data@MPeff, DLMdata@CV_Ind, DLMdata@Iref"
  ind <- max(1, (length(Data@Year) - yrsmth + 1)):length(Data@Year)
  deltaI <- mean(Data@Ind[x, ind])/Data@Iref[x]
  if (deltaI < (1 - mc)) 
    deltaI <- 1 - mc
  if (deltaI > (1 + mc)) 
    deltaI <- 1 + mc
  
  Effort <- Data@MPeff[x] * deltaI * trlnorm(reps, 1, Data@CV_Ind[x])
  if (reps == 1) 
    Effort <- Data@MPeff[x] * deltaI
  Effort[Effort < 0.01] <- 0.01  # for simulations in case Effort goes negative
  Effort <- mean(Effort)
  Allocate <- 1
  Spatial <- c(1, 1)
  Vuln <- rep(NA, 3)
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
}
class(ITe5) <- "Input"




#' This input control sets the minimum length of fish caught to a fraction of
#' the length that maximises the biomass, Lopt.
#' 
#' This aim of this simple MP is restrict the catch of small fish to rebuild
#' the stock biomass towards the optimal length, Lopt, expressed in terms of
#' the growth parameters Lopt=b/(M/k+b) (Hordyk et al. (2014)
#' 
#' 
#' @usage minlenLopt1(x, Data, reps = 100, buffer = 0.1)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param buffer Parameter controlling the fraction of Lopt to set the minimum
#' length of fish caught: minlen=Lopt*(0.7+buffer).
#' @return The length at first caprture, LFC, and length at full selectivity
#' @author HF Geromont
#' @references Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N., and J.
#' Prince. 2014. Some explorations of the life history ratios to describe
#' length composition, spawning-per-recruit, and the spawning potential ratio
#' ICES Journal of Marine Science, doi:10.1093/icesjms/fst235.
#' @export minlenLopt1
minlenLopt1 <- function(x, Data, reps = 100, buffer = 0.1) {
  
  # Minimum length MPs: Fix length-at-full-selectivity to 0.8*Lopt and
  # set length-at-first-capture 10% below LFs
  
  dependencies = "Data@MPeff, Data@vbLinf, Data@wlb, Data@Mort, Data@vbK"
  
  Lopt <- Data@vbLinf[x] * Data@wlb[x]/((Data@Mort[x]/Data@vbK[x]) + 
    Data@wlb[x])
  
  Allocate <- 1
  Spatial <- c(1, 1)
  newLFS <- Lopt * (0.7 + buffer)  # Lopt too precautionary, so set it to % below
  newLFC <- newLFS * 0.9
  Vuln <- c(newLFC, newLFS)
  Effort <- 1
  
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
  
}
class(minlenLopt1) <- "Input"




#' Effort MP: adjust effort up/down if mean length above/below Ltarget
#' 
#' Effort MP: adjust effort up/down if mean length above/below Ltarget
#' 
#' 
#' @usage EtargetLopt(x, Data, reps = 100, yrsmth=3, buffer=0.1)
#' @param x A position in data-limited methods data object
#' @param Data A data-limited methods data object
#' @param reps The number of TAC samples
#' @param yrsmth Number of years to calculate average length
#' @param buffer Parameter controlling the fraction of mean catch to set the
#' reference (or target) TAC level - acts as a precautionary buffer
#' @return An adjustment for fishing effort
#' @author HF Geromont
#' @export EtargetLopt
EtargetLopt <- function(x, Data, reps = 100, yrsmth = 3, buffer = 0.1) {
  
  # Effort MP: adjust effort up/down if mean length above/below Ltarget
  
  dependencies = "Data@Year, Data@ML, Data@L50, Data@MPeff, Data@vbLinf, Data@wlb, Data@Mort, Data@vbK"
  
  ind <- (length(Data@Year) - (yrsmth - 1)):length(Data@Year)  # recent 3 years
  Lrecent <- mean(Data@ML[ind])
  Lopt <- Data@vbLinf[x] * Data@wlb[x]/((Data@Mort[x]/Data@vbK[x]) + 
    Data@wlb[x])
  ratio <- Lrecent/Lopt
  
  Allocate <- 1
  Spatial <- c(1, 1)
  Vuln <- c(NA, 2)
  w <- 0.5
  Effort <- (1 - buffer) * (w + (1 - w) * ratio)
  
  out <- c(Allocate, Effort, Spatial, Vuln)
  out
  
}
class(EtargetLopt) <- "Input"





