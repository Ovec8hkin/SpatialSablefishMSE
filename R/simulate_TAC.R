#' Simulate Total Allowable Catch Projections
#' 
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
#' @param abc_tac_reduction proportional reduction from ABC to TAC
#' @param tac_land_reduction proportional reduction from TAC to landings
#'
#' @export simulate_TAC
#'
#' @example
#'
simulate_TAC <- function(hcr_F, naa, recruitment, joint_sel, dem_params, hist_tac, hcr_options, options){
    proj_faa <- joint_sel*hcr_F
    proj_N_new <- afscOM::simulate_population(naa, proj_faa, recruitment, dem_params, options=list())
    abc <- afscOM::baranov(hcr_F, proj_N_new$naa, dem_params$waa, dem_params$mort, joint_sel)
    tac <- abc * options$abc_tac_reduction

    # Implements symmetric stability constraints
    if(!is.na(hcr_options$max_stability) & !is.na(hist_tac)){
        max_tac <- hist_tac*(1+hcr_options$max_stability)
        min_tac <- hist_tac*(1-hcr_options$max_stability)
        if(tac > max_tac){
            tac <- max_tac
        }else if(tac < min_tac){
            tac <- min_tac
        }
    }

    # Implements a maximum tac cap
    if(!is.na(hcr_options$harvest_cap)){
        tac <- ifelse(tac > hcr_options$harvest_cap, hcr_options$harvest_cap, tac)
    }

    land <- tac * options$tac_land_reduction

    return(afscOM::listN(abc, tac, land))
}
