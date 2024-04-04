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
    #Try and improve the optimsation running the model for two additional Newton Raphson iterations
    try_improve = tryCatch(expr =
                            for(i in 1:2) {
                                g = as.numeric(my_model$gr(mle_optim$par))
                                h = optimHess(mle_optim$par, fn = my_model$fn, gr = my_model$gr)
                                mle_optim$par = mle_optim$par - solve(h,g)
                                mle_optim$objective = my_model$fn(mle_optim$par)
                            }
                        , error = function(e){e})
    
    mle_report = my_model$report(my_model$env$last.par.best)
    return(list(report=mle_report, model=my_model, opt=mle_optim))
}
