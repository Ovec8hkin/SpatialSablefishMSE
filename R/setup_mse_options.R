#' Setup MSE Options List
#' 
#' Set MSE options list items to reasonable default values.
#'
#' @export setup_mse_options
#'
#' @example
#'
setup_mp_options <- function(){

    return(
        list(
            hcr = NULL,
            ref_points = list(
                spr_target = 0.40
            ),
            management = list(
                abc_tac_reduction = 1,
                tac_land_reduciton = 1
            )
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
