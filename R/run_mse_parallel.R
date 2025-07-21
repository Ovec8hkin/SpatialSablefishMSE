#' Run Multiple MSE Simulations in Parallel
#' 
#' @param nsims number of simulations to run
#' @param seeds simulation seeds to use
#' @param om operating model object
#' @param hcr management procedure object
#' @param mse_options MSE options object
#' @param diagnostics whether to return EM model objects for diagnostic testing
#' 
#' @return list of MSE outputs (NAA, CAA, FAA, etc.)
#' 
#' @export run_mse_parallel
#' 
run_mse_parallel <- function(nsims, seeds, om, hcr, mse_options, diagnostics=FALSE, ...){

    nyears <- mse_options$n_proj_years+mse_options$n_spinup_years
    outputs <- setup_output_arrays(mse_options$n_proj_years+mse_options$n_spinup_years, nsims, seeds, mse_options$n_spinup_years)

    cores <- min(parallel::detectCores()-2, nsims)
    cl <- parallel::makeCluster(cores, outfile="")
    registerDoParallel(cl)

    out <- pbapply::pblapply(1:nsims, function(s, om, hcr, seeds, options){
        suppressMessages({
            library(tidyverse)
            library(TMB)
            library(devtools)
            library(abind)
            # library(SPoCK)
            library(afscOM)
            lapply(list.files("R", full.names = TRUE), source)
            # devtools::load_all("~/Desktop/Projects/afscOM")
        })
        
        seed <- seeds[s]
        mse <- run_mse(om=om, mp=hcr, mse_options=options, seed=seed)
        return(mse)

    }, om=om, hcr=hcr, seeds=seeds, options=mse_options, cl=cl)

    parallel::stopCluster(cl)

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
        outputs$exp_land[,,,,,s] <- mse$exp_land
        outputs$hcr_f[,,,,s] <- mse$hcr_f
        outputs$abc[,,,,,s] <- mse$abc
        outputs$tac[,,,,,s] <- mse$tac
        outputs$global_rec_devs[,,,,s] <- mse$global_rec_devs

        if(diagnostics){
            outputs$survey_obs$ll_rpn[,,,,s] <- mse$survey_obs$rpns[,,,,1]
            outputs$survey_obs$ll_rpw[,,,,s] <- mse$survey_obs$rpws[,,,,1]
            outputs$survey_obs$tw_rpw[,,,,s] <- mse$survey_obs$rpws[,,,,2]
            outputs$survey_obs$ll_acs[,,,,s] <- mse$survey_obs$acs[,,,,3]
            outputs$survey_obs$tw_acs[,,,,s] <- mse$survey_obs$acs[,,,,4]
            outputs$survey_obs$fxfish_acs[,,,,s] <- mse$survey_obs$acs[,,,,1]
            outputs$survey_obs$twfish_acs[,,,,s] <- mse$survey_obs$acs[,,,,2]

            if(mse_options$run_estimation){
                outputs$model_outs[(1:(nyears-mse_options$n_spinup_years+1))+((s-1)*(nyears-mse_options$n_spinup_years+1))] <- mse$model_outs[(1:(nyears-mse_options$n_spinup_years+1))]
            }
        }

    }

    outputs$dem_params <- om$dem_params
    outputs$mp <- mp

    return(outputs)

}

setup_output_arrays <- function(nyears, nsims, seeds, spinup_years){
    dimension_names <- list(
        "time" = 1:nyears,
        "age"  = 2:31,
        "sex"  = c("F", "M"),
        "region" = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
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
    faa_est     = array(NA, dim=c(nyears, nages, nsexes, 1, nfleets, nsims), dimnames=list("time"=1:(nyears), "age"=2:31, "sex"=c("F", "M"), "region"=c("Alaska"), "fleet"=c("Fixed", "Trawl"), "sim"=seeds))
    abc         = array(NA, dim=c(nyears+1, 1, 1, nregions, nfleets, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl"), "sim"=seeds))
    tac         = array(NA, dim=c(nyears+1, 1, 1, nregions, nfleets, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl"), "sim"=seeds))
    exp_land    = array(NA, dim=c(nyears+1, 1, 1, nregions, nfleets, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl"), "sim"=seeds))
    hcr_f       = array(NA, dim=c(nyears+1, 1, 1, 1, nsims), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska", "sim"=seeds))
    out_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sim"=seeds))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "sim"=seeds))
    naa_est     = array(NA, dim=c(nyears, nages, nsexes, 1, nsims), dimnames=list("time"=1:(nyears), "age"=2:31, "sex"=c("F", "M"), "region"=c("Alaska"), "sim"=seeds))
    global_rec_devs       = array(NA, dim=c(nyears-spinup_years, 1, 1, nregions, nsims), dimnames=list("time"=1:(nyears-spinup_years), 1, 1,"region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "sim"=seeds))

    survey_obs <- list(
        ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "sim"=seeds)),
        ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "sim"=seeds)),
        tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "sim"=seeds)),
        ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        tw_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")]),
        twfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=dimension_names[c("time", "age", "sex", "region", "sim")])
    )

    # model_outs = list(
    #     mods = array(list(), dim=c(nyears-spinup_years+1, nsims)),
    #     fits = array(list(), dim=c(nyears-spinup_years+1, nsims)),
    #     reps = array(list(), dim=c(nyears-spinup_years+1, nsims))
    # )
    model_outs <- vector("list", length=(nyears-spinup_years)*nsims)

    return(afscOM::listN(land_caa, disc_caa, caa, faa, faa_est, naa, naa_est, out_f, global_rec_devs, exp_land, abc, tac, hcr_f, survey_obs, model_outs))

}

