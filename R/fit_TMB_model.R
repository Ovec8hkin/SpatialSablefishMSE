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
fit_TMB_model <- function(data, parameters, fix_pars=NA){
    
    par_map <- list(
        # Turn off japanese fishery/survey and CPUE parameters
        ln_srv_jap_ll_sel_pars          = factor(rep(NA, length(parameters$ln_srv_jap_ll_sel_pars))),
        ln_srv_jap_fishery_ll_sel_pars  = factor(rep(NA, length(parameters$ln_srv_jap_fishery_ll_sel_pars))),
        ln_srv_jap_ll_q                 = factor(rep(NA, length(parameters$ln_srv_jap_ll_q))),
        ln_srv_jap_fishery_ll_q         = factor(rep(NA, length(parameters$ln_srv_jap_fishery_ll_q))),
        ln_ll_cpue_q                    = factor(rep(NA, length(parameters$ln_ll_cpue_q))),
    
        # Share selectivity parameters
        ln_trwl_sel_pars    = factor(c(1, 2, 3, 2)),
        ln_ll_sel_pars      = factor(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 4, 10, 11)),

        # Turn off mortality related parameters
        ln_M_year_devs  = factor(rep(NA, length(parameters$ln_M_year_devs))),
        ln_M_age_devs   = factor(rep(NA, length(parameters$ln_M_age_devs)))
    )

    # Need to do this for the selectivity estimation to not be weird
    data$ages <- as.double(1:30)

    if(!all(is.na(fix_pars))){
        parameters[names(fix_pars)] <- fix_pars
        par_map[names(fix_pars)] <- lapply(seq_along(fix_pars), \(x) factor(rep(NA, length(fix_pars[[x]]))))
    }
    
    dyn.load(dynlib("inst/CurrentAssessment"))
    my_model = TMB::MakeADFun(data = data,
                          parameters = parameters,
                          map = par_map,
                          DLL = "CurrentAssessment",
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
