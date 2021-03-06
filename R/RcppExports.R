# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Internal estimation function for LBSPR MP
#'
#' @param SL50 Length at 50 percent selectivity
#' @param SL95 Length at 95 percent selectivity
#' @param FM Ratio of apical fishing mortality to natural mortality
#' @param nage Number of pseudo age-classes
#' @param nlen Number of length bins
#' @param CVLinf CV of length-at-age
#' @param LenBins Vector of length bins
#' @param LenMids Vector of mid-points of length bins
#' @param MK Ratio of M/K
#' @param Linf Asymptotic length
#' @param rLens Vector of relative length at ate
#' @param Prob ALK 
#' @param Ml Maturity at age vector
#' @param L50 Length at 50 percent maturity
#' @param L95 Length at 95 percent maturity
#' @param Beta Exponent of the length-weight relationship
#' @author A. Hordyk
#' @useDynLib DLMtool
#' @keywords internal
#' @export
LBSPRgen <- function(SL50, SL95, FM, nage, nlen, CVLinf, LenBins, LenMids, MK, Linf, rLens, Prob, Ml, L50, L95, Beta) {
    .Call('_DLMtool_LBSPRgen', PACKAGE = 'DLMtool', SL50, SL95, FM, nage, nlen, CVLinf, LenBins, LenMids, MK, Linf, rLens, Prob, Ml, L50, L95, Beta)
}

LBSPRopt <- function(pars, CAL, nage, nlen, CVLinf, LenBins, LenMids, MK, Linf, rLens, Prob, Ml, L50, L95, Beta) {
    .Call('_DLMtool_LBSPRopt', PACKAGE = 'DLMtool', pars, CAL, nage, nlen, CVLinf, LenBins, LenMids, MK, Linf, rLens, Prob, Ml, L50, L95, Beta)
}

bhnoneq_LL <- function(stpar, year, Lbar, ss, Linf, K, Lc, nbreaks) {
    .Call('_DLMtool_bhnoneq_LL', PACKAGE = 'DLMtool', stpar, year, Lbar, ss, Linf, K, Lc, nbreaks)
}

