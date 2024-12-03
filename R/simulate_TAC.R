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
#' @param hist_tac TAC from the previous year
#' @param hcr_options list of additional management options such as stability
#' constraints and harvest caps
#' @param options list of additional ABC/TAC simulation options such as ABC-TAC
#' reduction levels and attainment functions.
#'
#' @export simulate_TAC
#'
#' @example
#'
simulate_TAC <- function(hcr_F, naa, recruitment, joint_sel, dem_params, hist_abc, hcr_options, options){
    proj_faa <- joint_sel*hcr_F
    proj_N_new <- afscOM::simulate_population(naa, proj_faa, recruitment, dem_params, options=list())
    abc <- afscOM::baranov(hcr_F, proj_N_new$naa, dem_params$waa, dem_params$mort, joint_sel)

    # Implements symmetric stability constraints
    if(!all(is.na(hcr_options$max_stability)) & !is.na(hist_abc)){
        if(length(hcr_options$max_stability) == 1){
            hcr_options$max_stability <- rep(hcr_options$max_stability, 2)
        }
        max_abc <- hist_abc*(1+hcr_options$max_stability[2])
        min_abc <- hist_abc*(1-hcr_options$max_stability[1])
        if(abc > max_abc){
            abc <- max_abc
        }else if(abc < min_abc){
            abc <- min_abc
        }
    }

    if(!is.list(options$abc_tac_reduction)){
        tac <- abc * options$abc_tac_reduction
    }else{
        tac <- abc*do.call(options$abc_tac_reduction$func, c(list(v=abc, naa=proj_N_new$naa), options$abc_tac_reduction$pars))
    }

    # Implements a maximum tac cap
    if(!is.na(hcr_options$harvest_cap)){
        tac <- ifelse(tac > hcr_options$harvest_cap, hcr_options$harvest_cap, tac)
    }

    # if(abc == 0 && tac > abc){
    #     tac <- 0
    # }

    if(!is.list(options$tac_land_reduction)){
        land <- tac * options$tac_land_reduction
    }else{
        land <- tac*do.call(options$tac_land_reduction$func, c(list(v=tac), options$tac_land_reduction$pars))
    }
    

    return(afscOM::listN(abc, tac, land, proj_N_new$naa))
}
