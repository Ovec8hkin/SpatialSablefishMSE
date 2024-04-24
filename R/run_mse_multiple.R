#' Run Multiple MSE Simulations
#' #'
#' Wrapper function around `run_mse` that handle running and
#' compiling the outputs of multiple MSE closed-loop simulations,
#' across different random seeds.
#'
#' @param nsims the number of unique simulations to run
#' @param seeds vector of random seeds to use (length(seeds) == nsims)
#' @param nyears number of years in each simulation
#' @param ... additional parameter sto pass to the `run_mse` call
#'
#' @return list of arrays of containing OM results from each simulation
#' @export run_mse_multiple
#' 
#' @example
#'
run_mse_multiple <- function(nsims, seeds, om, hcr, nyears, spinup_years=64, ...){

    dimension_names <- list(
        "time" = 1:nyears,
        "age"  = 2:31,
        "sex"  = c("F", "M"),
        "region" = "alaska",
        "fleet" = c("Fixed", "Trawl"),
        "sim" = seeds
    )

    nages <- length(dimension_names[["age"]])
    nsexes <- length(dimension_names[["sex"]])
    nregions <- length(dimension_names[["region"]])
    nfleets <- length(dimension_names[["fleet"]])

    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    faa_est     = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    abc         = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    tac         = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    exp_land    = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    hcr_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds))
    out_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=seeds))
    naa_est     = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=seeds))

    survey_obs <- list(
        ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds)),
        ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds)),
        tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds)),
        ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")])
    )

    model_outs = list(
        mods = array(list(), dim=c(nyears-spinup_years+1, nsims)),
        fits = array(list(), dim=c(nyears-spinup_years+1, nsims)),
        reps = array(list(), dim=c(nyears-spinup_years+1, nsims))
    )

    for(s in 1:nsims){
        seed <- seeds[s]
        mse <- run_mse(om, hcr, ..., nyears_input=nyears, spinup_years=spinup_years, seed=seed, file_suffix = seed)

        land_caa[,,,,,s] <- mse$land_caa
        disc_caa[,,,,,s] <- mse$disc_caa
        caa[,,,,,s] <- mse$caa
        faa[,,,,,s] <- mse$faa
        faa_est[,,,,,s] <- mse$faa_est
        naa[,,,,s] <- mse$naa
        naa_est[,,,,s] <- mse$naa_est
        out_f[,,,,s] <- mse$out_f
        exp_land[,,,,s] <- mse$exp_land
        hcr_f[,,,,s] <- mse$hcr_f
        abc[,,,,s] <- mse$abc
        tac[,,,,s] <- mse$tac

        survey_obs$ll_rpn[,,,,s] <- mse$survey_obs$ll_rpn
        survey_obs$ll_rpw[,,,,s] <- mse$survey_obs$ll_rpw
        survey_obs$tw_rpw[,,,,s] <- mse$survey_obs$tw_rpw
        survey_obs$ll_acs[,,,,s] <- mse$survey_obs$ll_acs
        survey_obs$fxfish_acs[,,,,s] <- mse$survey_obs$fxfish_acs

        model_outs$mods[,s] <- mse$model_outs$mods
        model_outs$fits[,s] <- mse$model_outs$fits
        model_outs$reps[,s] <- mse$model_outs$reps

    }

    return(afscOM::listN(land_caa, disc_caa, caa, faa, faa_est, naa, naa_est, out_f, exp_land, abc, tac, hcr_f, survey_obs, model_outs))

}
