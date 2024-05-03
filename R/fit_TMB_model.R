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
fit_TMB_model <- function(data, parameters, model_name="CurrentAssessment", do_newton_steps=FALSE, fix_pars=NA, recompile_model=FALSE){
    
    if(!(model_name %in% c("CurrentAssessment", "CurrentAssessmentDisaggregated"))){
        stop("Invalid model_name")
    }

    par_map <- list(
        # Turn off japanese fishery/survey and CPUE parameters
        ln_srv_jap_ll_sel_pars          = factor(rep(NA, length(parameters$ln_srv_jap_ll_sel_pars))),
        ln_srv_jap_fishery_ll_sel_pars  = factor(rep(NA, length(parameters$ln_srv_jap_fishery_ll_sel_pars))),
        ln_srv_jap_ll_q                 = factor(rep(NA, length(parameters$ln_srv_jap_ll_q))),
        ln_srv_jap_fishery_ll_q         = factor(rep(NA, length(parameters$ln_srv_jap_fishery_ll_q))),
        ln_ll_cpue_q                    = factor(rep(NA, length(parameters$ln_ll_cpue_q))),

        # Selectivity parameter sharing
        # ln_trwl_sel_pars    = factor(c(1, 2, 3, 2)),

        # Turn off mortality related parameters
        ln_M_year_devs  = factor(rep(NA, length(parameters$ln_M_year_devs))),
        ln_M_age_devs   = factor(rep(NA, length(parameters$ln_M_age_devs)))
    )

    # Need to do this for the selectivity estimation to not be weird
    data$ages <- as.double(1:30)

    if(!all(is.na(fix_pars))){
        parameters[names(fix_pars)] <- lapply(seq_along(fix_pars), \(x) as.vector(fix_pars[[x]]))
        par_map[names(fix_pars)] <- lapply(seq_along(fix_pars), \(x) factor(rep(NA, length(fix_pars[[x]]))))
    }

    if(recompile_model){
        file.remove(paste0("inst/", model_name, ".o"))
        file.remove(paste0("inst/", model_name, ".so"))
        TMB::compile(paste0("inst/", model_name, ".cpp"))
    }
    
    dyn.load(dynlib(paste0("inst/", model_name)))
    my_model = TMB::MakeADFun(data = data,
                          parameters = parameters,
                          map = par_map,
                          DLL = model_name,
                          silent = TRUE)

    mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))
    
    #Try and improve the optimsation running the model for two additional Newton Raphson iterations
    if(do_newton_steps){
        tryCatch(expr =
            for(i in 1:2) {
                g = as.numeric(my_model$gr(mle_optim$par))
                h = optimHess(mle_optim$par, fn = my_model$fn, gr = my_model$gr)
                mle_optim$par = mle_optim$par - solve(h,g)
                mle_optim$objective = my_model$fn(mle_optim$par)
            }
        , error = function(e){e})
    }
    
    mle_report = my_model$report(my_model$env$last.par.best)
    return(list(report=mle_report, model=my_model, opt=mle_optim))
}
