#' Compute F from Tier 3 HCR
#'
#' @param ssb Value for SSB
#' @param B_ref Value for biomass based reference point (e.g., B40)
#' @param F_ref Value for fishing mortality based reference point (e.g., F40)
#'
#' @return fishing mortality rate
#' @export npfmc_tier3_F
#'
#' @examples
npfmc_tier3_F <- function(ssb, B_ref, F_ref){
    x <- ssb/B_ref
    return(
        threshold_f(x, f_min=0, f_max=F_ref, lrp=0.05, u=1)
    )
}

#' Compute F from a constant F rule
#'
#' @param F fishing mortality rate to use
#'
#' @return fishing mortality rate
#' @export constant_F
#'
#' @example
#'
constant_F <- function(F){
    return(F)
}

#' Compute F following a threshold rule, with a cap at B/B_ref = 1
#'
#' @param ssb Value for SSB
#' @param B_ref Value for biomass based reference point (e.g., B40)
#' @param F_ref Value for fishing mortality based reference point (e.g., F40)
#' @param alpha Cutoff to set F at 0
#'
#' @return fishing mortality rate
#' @export threshold_cap
#'
#' @examples
#' 
threshold_cap <- function(ssb,B_ref,F_ref,alpha) {
  if(ssb/B_ref >= 1)  F <- F_ref  / (ssb/B_ref)  # if stock status >= 1
  if(ssb/B_ref > alpha && ssb/B_ref < 1) F <- F_ref * ((ssb/B_ref)-alpha)/(1-alpha) # if stock status > alpha & stock status < 1
  if(ssb/B_ref <= alpha) F <- 0 # if stock stats <= alpha
  return(F)
}

#' Generic threshold F harvest control rule
#' #'
#' A generic 4 parameters threshold F HCR. 
#'
#' @param x input to HCR (e.g. SSB or, SSB/B40)
#' @param f_min F to use when x < lrp
#' @param f_max F to use when x > urp
#' @param lrp value of x below which f_min applies
#' @param urp value of x above which f_max applies
#'
#' @return fishing mortality rate
#' @export threshold_f
#'
#' @example
#'
threshold_f <- function(x, f_min, f_max, lrp, urp){
    if(x >= urp)  F <- f_max # if stock status >= 1
    if(x > lrp && x < urp) F <- F_ref * ((x-lrp)/(urp-lrp) # if stock status > alpha & stock status < 1
    if(x <= lrp) F <- f_min # if stock stats <= alpha
    return(F)
}