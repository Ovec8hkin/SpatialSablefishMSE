#' Setup Management Procedure Object
#' 
#' Set MP object list items to reasonable default values.
#'
#' @export setup_mse_options
#'
#' @example
#'
setup_mp_options <- function(){

    return(
        list(
            hcr = setup_hcr_options(),
            ref_points = list(
                spr_target = 0.40
            ),
            management = list(
                abc_tac_reduction = 1,
                tac_land_reduciton = 1
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
            n_proj_years = 100,
            n_spinup_years = 64,
            recruitment_start_year = 64,
            run_estimation = TRUE
        )
    )
}
