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
    
    ## 1. Use harvest control rule to determine ABC including
    ## all stability constraints and harvest caps
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

    # Implements a maximum tac cap
    if(!is.na(hcr_options$harvest_cap)){
        abc <- ifelse(abc > hcr_options$harvest_cap, hcr_options$harvest_cap, abc)
    }

    ## 2. Apportion ABC to different regions using supplied 
    ## apportionment scheme
    abc_reg <- abc*options$abc_region_apportionment

    ## 3. Allocate regional ABCs to fleets within regions based
    ## on supplied region-fleet splits
    abc_regflt <- abc_reg*options$abc_regflt_apportionment

    ## 4. Allow for TAC to differ from ABC based on supplied
    ## values.
    tac_regflt <- tac <- abc_regflt * options$abc_tac_regflt_reduction
    
    # if(!is.list(options$abc_tac_reduction)){
    #     tac <- abc_regflt * options$abc_tac_regflt_reduction
    # }else{
    #     tac <- abc*do.call(options$abc_tac_reduction$func, c(list(v=abc, naa=proj_N_new$naa), options$abc_tac_reduction$pars))
    # }

    ## 6. Allow for region-fleet specific TAC utilizations based
    ## on supplied values.
    land_regflt <- tac_regflt * options$regflt_tac_utilization

    # if(!is.list(options$tac_land_reduction)){
    #     land <- tac * options$tac_land_reduction
    # }else{
    #     land <- tac*do.call(options$tac_land_reduction$func, c(list(v=tac), options$tac_land_reduction$pars))
    # }
    
    # abc, tac, land are now all of dimensions [1, nages, nsexes, nregions, nfleets]

    # return(afscOM::listN(abc, tac, land, proj_N_new$naa))
    return(list(
        abc = abc_regflt,
        tac = tac_regflt,
        land = land_regflt,
        proj_N_new = proj_N_new$naa
    ))
}
