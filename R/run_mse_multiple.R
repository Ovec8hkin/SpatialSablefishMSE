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
#' @export list of arrays of containing OM results from each simulation
#'
#' @example
#'
run_mse_multiple <- function(nsims, seeds, nyears, ...){

    dimension_names <- list(
        "time" = 1:nyears,
        "age"  = 2:31,
        "sex"  = c("F", "M"),
        "region" = "alaska",
        "fleet" = c("Fixed", "Trawl"),
        "sim" = seeds
    )

    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    tac         = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds))
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

    for(s in 1:nsims){
        seed <- seeds[s]
        mse <- run_mse(om=om, hcr=hcr, nyears_input=nyears, seed=seed)
    
        land_caa[,,,,,s] <- mse$land_caa
        disc_caa[,,,,,s] <- mse$disc_caa
        caa[,,,,,s] <- mse$caa
        faa[,,,,,s] <- mse$faa
        naa[,,,,s] <- mse$naa
        naa_est[,,,,s] <- mse$naa_est
        out_f[,,,,s] <- mse$out_f
        tac[,,,,s] <- mse$tac
        hcr_f[,,,,s] <- mse$hcr_f

        survey_obs$ll_rpn[,,,,s] <- mse$survey_obs$ll_rpn
        survey_obs$ll_rpw[,,,,s] <- mse$survey_obs$ll_rpw
        survey_obs$tw_rpw[,,,,s] <- mse$survey_obs$tw_rpw
        survey_obs$ll_acs[,,,,s] <- mse$survey_obs$ll_acs
        survey_obs$fxfish_acs[,,,,s] <- mse$survey_obs$fxfish_acs

    }

    return(afscOM::listN(land_caa, disc_caa, caa, faa, naa, naa_est, out_f, tac, hcr_f, survey_obs))

}
