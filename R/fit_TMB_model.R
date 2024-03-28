#' Fit the SpatialSablefishAssessment TMB Model
#' #'
#' Use Template Model Builder (TMB) to fit the assessmnte model
#' from `SpatialSablefishAssessment`. Note that the 'ln_M_year_devs'
#' and 'ln_M_age_devs' parameters have been turned off.
#'
#' @param data data list to be provided to TMB (from `format_em_data(...)`)
#' @param parameters parameter list to be provided to TMB (from `format_em_data(...)`)
#'
#' @export fit_TMB_model
#'
#' @example
#'
fit_TMB_model <- function(data, parameters){
    par_map <- list(
        ln_M_year_devs = factor(rep(NA, length(parameters$ln_M_year_devs))),
        ln_M_age_devs = factor(rep(NA, length(parameters$ln_M_age_devs)))
    )
    my_model = TMB::MakeADFun(data = data,
                          parameters = parameters,
                          map = par_map,
                          DLL = "SpatialSablefishAssessment_TMBExports",
                          silent = TRUE)

    mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))
    mle_report = my_model$report(mle_optim$par)
    return(mle_report)
}
