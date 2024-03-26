fit_TMB_model <- function(data, parameters){
    my_model = TMB::MakeADFun(data = data,
                          parameters = parameters,
                          DLL = "SpatialSablefishAssessment_TMBExports",
                          silent = TRUE)

    mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))
    mle_report = my_model$report(mle_optim$par)
    return(mle_report)
}
