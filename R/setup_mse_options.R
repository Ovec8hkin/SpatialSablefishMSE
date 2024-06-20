#' Setup MSE Options List
#' 
#' Set MSE options list items to reasonable default values.
#'
#' @export setup_mse_options
#'
#' @example
#'
setup_mse_options <- function(){

    return(
        list(
            hcr = NULL,
            ref_points = list(
                spr_target = 0.40
            ),
            management = list(
                abc_tac_reduction = 1,
                tac_land_reduction = 1
            )
        )
    )

}
