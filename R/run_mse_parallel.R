run_mse_parallel <- function(nsims, seeds, om, hcr, mse_options, nyears, diagnostics=FALSE, ...){

    outputs <- setup_output_arrays(nyears, nsims, seeds, mse_options$n_spinup_years)

    cores <- min(parallel::detectCores()-2, nsims)
    cl <- parallel::makeCluster(cores, outfile="")
    registerDoParallel(cl)

    out <- pbapply::pblapply(1:nsims, function(s, om, hcr, nyears, seeds, options){
        suppressMessages({
            library(tidyverse)
            library(TMB)
            library(devtools)
            library(abind)

            lapply(list.files("R", full.names = TRUE), source)
            devtools::load_all("~/Desktop/Projects/afscOM")
        })
        
        seed <- seeds[s]
        mse <- run_mse(om=om, mp=hcr, mse_options=options, nyears_input=nyears, seed=seed, file_suffix = seed)
        return(mse)

    }, om=om, hcr=hcr, nyears=nyears, seeds=seeds, options=mse_options, cl=cl)

    stopCluster(cl)

    for(s in 1:nsims){

        mse <- out[[s]]

        outputs$land_caa[,,,,,s] <- mse$land_caa
        outputs$disc_caa[,,,,,s] <- mse$disc_caa
        outputs$caa[,,,,,s] <- mse$caa
        outputs$faa[,,,,,s] <- mse$faa
        outputs$faa_est[,,,,,s] <- mse$faa_est
        outputs$naa[,,,,s] <- mse$naa
        outputs$naa_est[,,,,s] <- mse$naa_est
        outputs$out_f[,,,,s] <- mse$out_f
        outputs$exp_land[,,,,s] <- mse$exp_land
        outputs$hcr_f[,,,,s] <- mse$hcr_f
        outputs$abc[,,,,s] <- mse$abc
        outputs$tac[,,,,s] <- mse$tac

        if(diagnostics){
            outputs$survey_obs$ll_rpn[,,,,s] <- mse$survey_obs$rpns[,,,,1]
            outputs$survey_obs$ll_rpw[,,,,s] <- mse$survey_obs$rpws[,,,,1]
            outputs$survey_obs$tw_rpw[,,,,s] <- mse$survey_obs$rpws[,,,,2]
            outputs$survey_obs$ll_acs[,,,,s] <- mse$survey_obs$acs[,,,,3]
            outputs$survey_obs$tw_acs[,,,,s] <- mse$survey_obs$acs[,,,,4]
            outputs$survey_obs$fxfish_acs[,,,,s] <- mse$survey_obs$acs[,,,,1]
            outputs$survey_obs$twfish_acs[,,,,s] <- mse$survey_obs$acs[,,,,2]

            if(om$model_options$run_estimation){
                outputs$model_outs$mods[,s] <- mse$model_outs$mods
                outputs$model_outs$fits[,s] <- mse$model_outs$fits
                outputs$model_outs$reps[,s] <- mse$model_outs$reps
            }
        }

    }

    return(outputs)

}

setup_output_arrays <- function(nyears, nsims, seeds, spinup_years){
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
    faa_est         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    abc         = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    tac         = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    exp_land    = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    hcr_f       = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    out_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=seeds))
    naa_est     = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=seeds))

    survey_obs <- list(
        ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds)),
        ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds)),
        tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds)),
        ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        tw_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        twfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")])
    )

    model_outs = list(
        mods = array(list(), dim=c(nyears-spinup_years+1, nsims)),
        fits = array(list(), dim=c(nyears-spinup_years+1, nsims)),
        reps = array(list(), dim=c(nyears-spinup_years+1, nsims))
    )

    return(afscOM::listN(land_caa, disc_caa, caa, faa, faa_est, naa, naa_est, out_f, exp_land, abc, tac, hcr_f, survey_obs, model_outs))

}