#' Simulate Total Allowable Catch Projections
#' #'
#' Determine the total allowable catch (TAC) in the next simulation
#' year based on the fishing mortality rate identified by a harvest
#' control rule (HCR). TAC is computed by projecting the population
#' forward one year under the specific level of F and computing the
#' the amount of catch that would be taken from that projected 
#' population structure based on the joint fishery selectivity curve.
#'
#' @param hcr_F fishing mortality rate identified by HCR
#' @param naa numbers-at-age vector (dimensions [1, nages, nsexes, nregions])
#' @param recruitment projected recruitment in next year
#' @param join_sel joint fishery selectivity
#' @param dem_params demographic parameter matrices subsetted to 1 year
#'
#' @export simulate_TAC
#'
#' @example
#'
simulate_TAC <- function(hcr_F, naa, recruitment, joint_sel, dem_params){
    proj_faa <- joint_sel*hcr_F
    proj_N_new <- afscOM::simulate_population(naa, proj_faa, recruitment, dem_params, options=list())
    tac <- afscOM::baranov(hcr_F, proj_N_new$naa, dem_params$waa, dem_params$mort, joint_sel)

    return(tac)
}
