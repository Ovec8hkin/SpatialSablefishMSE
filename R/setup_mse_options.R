#' Setup Management Procedure Object
#' 
#' Set MP object list items to reasonable default values.
#'
#' @export setup_mp_options
#'
#' @example
#'
setup_mp_options <- function(){

    return(
        list(
            hcr = setup_hcr_options(),
            ref_points = list(
                spr_target = c(0.40, 0.40)
            ),
            management = list(
                abc_regflt_apportionment = abc_regionfleet_allocation,
                abc_tac_regflt_reduction = array(1, dim=c(5, 2), dimnames=dimnames(abc_regionfleet_allocation)),
                regflt_tac_utilization = average_regflt_tac_utilization
            ),
            apportionment = list(
                func = rpw_moving_average,
                pars = list(window_size=5)
            ),       
            survey_frequency = 1,
            assessment_frequency = 1
        )
    )

}

#' Setup HCR Objects
#' 
#' Set HCR object list items to NULL defaults
#' 
#' @export setup_hcr_options
#' 
setup_hcr_options <- function(){
    return(
        list(
            func = NULL,
            extra_pars = NA,
            extra_options = list(
                max_stability = NA,
                harvest_cap = NA
            ),
            units = NULL
        )
    )
}

setup_mse_options <- function(){
    return(
        list(
            n_proj_years = 75,
            n_spinup_years = 54,
            recruitment_start_year = 54,
            run_estimation = TRUE
        )
    )
}

#' ABC Alloaction Rate to Fleets within Regions
#' 
#' Values represent the proportion of the regional ABC that is allocated to each fleet.
#' Values come from the GOA Groudfish FMP and BSAI Groundfish FMP.
#' 
#' @export abc_regionfleet_allocation
#' 
abc_regionfleet_allocation <- array(
    c(0.5, 0.75, 0.8, 0.8, 0.95, 
      0.5, 0.25, 0.2, 0.2, 0.05), # Fixed and Trawl fleets
    dim = c(5, 2),
    dimnames = list(
        "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"),
        "fleet"=c("Fixed", "Trawl")
    )
)
