#' Title Compute F from Tier 3 HCR
#'
#' @param ssb Value for SSB
#' @param B_ref Value for biomass based reference point (e.g., B40)
#' @param F_ref Value for fishing mortality based reference point (e.g., F40)
#' @param alpha Cutoff to set F at 0
#'
#' @return
#' @export
#'
#' @examples
npfmc_tier3_F <- function(ssb, B_ref, F_ref, alpha){
    if(ssb/B_ref >= 1)  F <- F_ref # if stock status >= 1
    if(ssb/B_ref > alpha && ssb/B_ref < 1) F <- F_ref * ((ssb/B_ref)-alpha)/(1-alpha) # if stock status > alpha & stock status < 1
    if(ssb/B_ref <= alpha) F <- 0 # if stock stats <= alpha
    return(F)
}

#' Title Compute F following a threshold rule, with a cap at B/B_ref = 1
#'
#' @param ssb Value for SSB
#' @param B_ref Value for biomass based reference point (e.g., B40)
#' @param F_ref Value for fishing mortality based reference point (e.g., F40)
#' @param alpha Cutoff to set F at 0
#'
#' @return
#' @export
#'
#' @examples
threshold_cap <- function(ssb,B_ref,F_ref,alpha) {
  if(ssb/B_ref >= 1)  F <- F_ref  / (ssb/B_ref)  # if stock status >= 1
  if(ssb/B_ref > alpha && ssb/B_ref < 1) F <- F_ref * ((ssb/B_ref)-alpha)/(1-alpha) # if stock status > alpha & stock status < 1
  if(ssb/B_ref <= alpha) F <- 0 # if stock stats <= alpha
  return(F)
}
