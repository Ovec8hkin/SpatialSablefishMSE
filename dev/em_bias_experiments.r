rm(list=ls())
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)
library(pbapply)

# library(afscOM) # may work but not certain

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM" 
devtools::load_all(afscOM_dir)

# spatialSablefish_dir <- "~/Desktop/Projects/SpatialSablefishAssessment" 
# devtools::load_all(spatialSablefish_dir)

source("R/reference_points.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/age_structure_stats.R")
source("R/data_utils.R")
source("R/data_processing.R")
source("R/run_mse.R")
source("R/run_mse_multiple.R")
source("R/format_em_data.R")
source("R/fit_TMB_model.R")

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons
sable_om <- readRDS("data/sablefish_om.RDS") # Read this saved OM from a file
sable_om$model_options$simulate_observations <- TRUE # Enable simulating observations
# Generate age comp samples as integers rather than proportions
# (this is necesarry for the SpatialSablefishAssessment TMB model)
sable_om$model_options$obs_pars$surv_ll$as_integers = TRUE
sable_om$model_options$obs_pars$surv_tw$as_integers = TRUE
sable_om$model_options$obs_pars$fish_fx$as_integers = TRUE
# Turn off estimation model
sable_om$model_options$run_estimation = FALSE

# OM with sex-spexific selectivity and small amount of observation error
om_sexed_small <- sable_om
om_sexed_small$model_options$obs_pars$surv_ll$ac_samps <- 150
om_sexed_small$model_options$obs_pars$surv_tw$ac_samps <- 150
om_sexed_small$model_options$obs_pars$fish_fx$ac_samps <- 150
om_sexed_small$model_options$obs_pars$surv_ll$rpn_cv <- 0.05
om_sexed_small$model_options$obs_pars$surv_ll$rpw_cv <- 0.05
om_sexed_small$model_options$obs_pars$surv_tw$rpw_cv <- 0.05

# OM with sex-aggregated selectivity and small amount of observation error
om_unsexed_small <- sable_om
om_unsexed_small$model_options$obs_pars$surv_ll$ac_samps <- 150
om_unsexed_small$model_options$obs_pars$surv_tw$ac_samps <- 150
om_unsexed_small$model_options$obs_pars$fish_fx$ac_samps <- 150
om_unsexed_small$model_options$obs_pars$surv_ll$rpn_cv <- 0.05
om_unsexed_small$model_options$obs_pars$surv_ll$rpw_cv <- 0.05
om_unsexed_small$model_options$obs_pars$surv_tw$rpw_cv <- 0.05
om_unsexed_small$dem_params$sel[,,2,,] <- om_unsexed_small$dem_params$sel[,,1,,]
om_unsexed_small$dem_params$surv_sel[,,2,,] <- om_unsexed_small$dem_params$surv_sel[,,1,,]

# OM with sex-spexific selectivity and small amount of observation error
om_sexed_medium <- sable_om
om_sexed_medium$model_options$obs_pars$surv_ll$ac_samps <- 100
om_sexed_medium$model_options$obs_pars$surv_tw$ac_samps <- 100
om_sexed_medium$model_options$obs_pars$fish_fx$ac_samps <- 100
om_sexed_medium$model_options$obs_pars$surv_ll$rpn_cv <- 0.1
om_sexed_medium$model_options$obs_pars$surv_ll$rpw_cv <- 0.1
om_sexed_medium$model_options$obs_pars$surv_tw$rpw_cv <- 0.1

# OM with sex-aggregated selectivity and medium amount of observation error
om_unsexed_medium <- sable_om
om_unsexed_medium$model_options$obs_pars$surv_ll$ac_samps <- 100
om_unsexed_medium$model_options$obs_pars$surv_tw$ac_samps <- 100
om_unsexed_medium$model_options$obs_pars$fish_fx$ac_samps <- 100
om_unsexed_medium$model_options$obs_pars$surv_ll$rpn_cv <- 0.1
om_unsexed_medium$model_options$obs_pars$surv_ll$rpw_cv <- 0.1
om_unsexed_medium$model_options$obs_pars$surv_tw$rpw_cv <- 0.1
om_unsexed_medium$dem_params$sel[,,2,,] <- om_unsexed_medium$dem_params$sel[,,1,,]
om_unsexed_medium$dem_params$surv_sel[,,2,,] <- om_unsexed_medium$dem_params$surv_sel[,,1,,]

# OM with sex-spexific selectivity and large amount of observation error
om_sexed_large <- sable_om
om_sexed_large$model_options$obs_pars$surv_ll$ac_samps <- 50
om_sexed_large$model_options$obs_pars$surv_tw$ac_samps <- 50
om_sexed_large$model_options$obs_pars$fish_fx$ac_samps <- 50
om_sexed_large$model_options$obs_pars$surv_ll$rpn_cv <- 0.20
om_sexed_large$model_options$obs_pars$surv_ll$rpw_cv <- 0.10
om_sexed_large$model_options$obs_pars$surv_tw$rpw_cv <- 0.10

# OM with sex-aggregated selectivity and large amount of observation error
om_unsexed_large <- sable_om
om_unsexed_large$model_options$obs_pars$surv_ll$ac_samps <- 50
om_unsexed_large$model_options$obs_pars$surv_tw$ac_samps <- 50
om_unsexed_large$model_options$obs_pars$fish_fx$ac_samps <- 50
om_unsexed_large$model_options$obs_pars$surv_ll$rpn_cv <- 0.20
om_unsexed_large$model_options$obs_pars$surv_ll$rpw_cv <- 0.10
om_unsexed_large$model_options$obs_pars$surv_tw$rpw_cv <- 0.10
om_unsexed_large$dem_params$sel[,,2,,] <- om_unsexed_large$dem_params$sel[,,1,,]
om_unsexed_large$dem_params$surv_sel[,,2,,] <- om_unsexed_large$dem_params$surv_sel[,,1,,]





#' 2. Define a harvest control rule (HCR) function to use to project TAC
#' in future years. Such function must take accept the following parameters:
#' ref_pts: a list of reference points as generated by `calculate_ref_points`,
#'          (B40, F35, and F40)
#' naa: a vector of numbers-at-age for a given year (should be [1, nages, nsexes, nregions])
#' dem_params: a list of demographic parameter values subsetted to a single year
#' 
#' HCR functions can accept as many additional parameters as required.
#' 
#' Below is an example implementation of the NPFMCs Tier 3a HCR     
tier3 <- function(ref_pts, naa, dem_params){
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
    )
}

#' 3. Run a lot of OM projections
#' A single OM simulation can be run using the `run_mse(...)`
#' function with om$model_options$run_estimation = FALSE, 
#' while multiple MSE simulations can be run (serially)
#' using the `run_mse_multiple(...)` function.
#' 
#' It is recommended to always use `run_mse_multiple(...)` even
#' when only a single MSE simulation is required.
set.seed(2023)
nsims <- 50
seeds <- sample(1:(1000*nsims), nsims)  # Draw 50 random seeds


mse_sexed_small     <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=160, om=om_sexed_small, hcr=tier3)
mse_unsexed_small   <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=160, om=om_unsexed_small, hcr=tier3)
mse_sexed_medium    <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=160, om=om_sexed_medium, hcr=tier3)
mse_unsexed_medium  <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=160, om=om_unsexed_medium, hcr=tier3)
mse_sexed_large     <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=160, om=om_sexed_large, hcr=tier3)
mse_unsexed_large   <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=160, om=om_unsexed_large, hcr=tier3)

####

#' 4. Apply the TMB fitting procedure to the observation data
#' Use fit_TMB_model(...) to fit the EM to the generated 
#' observations from each of the MSEs in the terminal year of
#' the OM projection.
#' 
#' This is being done in parallele using pbapply::pblapply(...) 
cores <- parallel::detectCores()-2
cl <- parallel::makeCluster(cores, outfile="")
registerDoParallel(cl)

assess_ssb_agg_sxsm <- pbapply::pblapply(1:50, function(s, sable_om, mse_tier3, om){
    assess_ssb <- run_EM(s, sable_om, mse_tier3, om)
    return(assess_ssb)
}, sable_om=om_sexed_small, mse_tier3=mse_sexed_small, om="sxsm", cl=cl)
assess_ssb_agg_sxsm <- bind_rows(assess_ssb_agg_sxsm)

assess_ssb_agg_usxsm <- pbapply::pblapply(1:50, function(s, sable_om, mse_tier3, om){
    assess_ssb <- run_EM(s, sable_om, mse_tier3, om)
    return(assess_ssb)
}, sable_om=om_unsexed_small, mse_tier3=mse_unsexed_small, om="usxsm", cl=cl)
assess_ssb_agg_usxsm <- bind_rows(assess_ssb_agg_usxsm)

assess_ssb_agg_sxlg <- pbapply::pblapply(1:50, function(s, sable_om, mse_tier3, om){
    assess_ssb <- run_EM(s, sable_om, mse_tier3, om)
    return(assess_ssb)
}, sable_om=om_sexed_large, mse_tier3=mse_sexed_large, om="sxlg", cl=cl)
assess_ssb_agg_sxlg <- bind_rows(assess_ssb_agg_sxlg)

assess_ssb_agg_usxlg <- pbapply::pblapply(1:50, function(s, sable_om, mse_tier3, om){
    assess_ssb <- run_EM(s, sable_om, mse_tier3, om)
    return(assess_ssb)
}, sable_om=om_unsexed_large, mse_tier3=mse_unsexed_large, om="usxlg", cl=cl)
assess_ssb_agg_usxlg <- bind_rows(assess_ssb_agg_usxlg)

assess_ssb_agg_sxmd <- pbapply::pblapply(1:50, function(s, sable_om, mse_tier3, om){
    assess_ssb <- run_EM(s, sable_om, mse_tier3, om)
    return(assess_ssb)
}, sable_om=om_sexed_medium, mse_tier3=mse_sexed_medium, om="sxmd", cl=cl)
assess_ssb_agg_sxmd <- bind_rows(assess_ssb_agg_sxmd)

assess_ssb_agg_usxmd <- pbapply::pblapply(1:50, function(s, sable_om, mse_tier3, om){
    assess_ssb <- run_EM(s, sable_om, mse_tier3, om)
    return(assess_ssb)
}, sable_om=om_unsexed_medium, mse_tier3=mse_unsexed_medium, om="usxmd", cl=cl)
assess_ssb_agg_usxmd <- bind_rows(assess_ssb_agg_usxmd)

stopCluster(cl)

#' 5. Process EM fits and compute relative error
assess_ssb_agg <- bind_rows(assess_ssb_agg_sxsm, assess_ssb_agg_sxlg, assess_ssb_agg_usxlg, assess_ssb_agg_usxmd, assess_ssb_agg_sxmd)

assess_ssb_agg %>% as_tibble() %>%
    mutate(rel_error = (SSB-om_ssb)/om_ssb) %>%
    group_by(Year, om) %>%
    median_qi(SSB, om_ssb, .width=c(0.50, 0.95)) %>%

    ggplot(aes(x=Year))+
        geom_lineribbon(aes(y=om_ssb, ymin=om_ssb.lower, ymax=om_ssb.upper))+
        scale_fill_brewer(palette="Blues")+
        geom_pointrange(aes(y=SSB, ymin=SSB.lower, ymax=SSB.upper), alpha=0.2, size=0., color="red")+
        scale_y_continuous(limits=c(0, 300))+
        facet_wrap(~om)+
        theme_bw() 

mu_error <- assess_ssb_agg %>% as_tibble() %>%
    filter(Year > 2022, Year < 2022+100) %>%
    mutate(rel_error = (SSB-om_ssb)/om_ssb) %>%
    group_by(om) %>%
    summarise(mu=median(rel_error))

assess_ssb_agg %>% as_tibble() %>%
    mutate(rel_error = (SSB-om_ssb)/om_ssb) %>%
    ungroup() %>%
    group_by(Year, om) %>%
    #mutate(mu=mean(rel_error)) %>%
    median_qi(rel_error, .width=c(0.50, 0.95)) %>%

    ggplot(aes(x=Year))+
        geom_lineribbon(aes(y=rel_error, ymin=.lower, ymax=.upper))+
        geom_hline(yintercept = 0)+
        geom_hline(data=mu_error, aes(yintercept = mu), color="red")+
        geom_vline(xintercept = 2023)+
        geom_text(data=mu_error, aes(x=2075, y=0.22, label=round(mu, 3)))+
        scale_fill_brewer(palette="Blues")+
        facet_wrap(~om)+
        #geom_pointrange(aes(y=SSB, ymin=SSB.lower, ymax=SSB.upper), color="red")+
        #scale_y_continuous(limits=c(0, 300))+
        theme_bw() 


#' Run Estimation Method 
#' 
#' Given an OM object and a succesfull MSE run, apply
#' the EM to the observation data generated by the MSE.
#'
#' @param s simulation seed
#' @param sable_om om object
#' @param mse_tier3 mse object
#' @param om string identifed for OM
#'
#' @export
#'
#' @example
#'
run_EM <- function(s, sable_om, mse_tier3, om){
    afscOM_dir <- "~/Desktop/Projects/afscOM"
    
    library(devtools) 
    devtools::load_all(afscOM_dir)

    source("R/format_em_data.R")
    source("R/fit_TMB_model.R")

    i=160-63
    y <- 63+i

    assess_inputs <- format_em_data(
                        nyears = y,
                        dem_params = afscOM::subset_dem_params(sable_om$dem_params, 64:y, d=1, drop=FALSE),
                        land_caa = afscOM::subset_matrix(mse_tier3$land_caa[64:y,,,,,s,drop=FALSE], 1, d=6, drop=TRUE),
                        survey_indices = afscOM::subset_dem_params(afscOM::subset_dem_params(mse_tier3$survey_obs, 64:y, d=1, drop=FALSE), s, d=5, drop=TRUE),
                        fxfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$fxfish_acs[63:(y-1),,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                        ll_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$ll_acs[63:(y-1),,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                        model_options = sable_om$model_options,
                        added_years = i,
                        file_suffix = s
                    )
    file.remove(paste0("data/sablefish_em_data_curr_",s,".RDS"))
    file.remove(paste0("data/sablefish_em_par_curr_",s,".RDS"))
    mod_out <- fit_TMB_model(assess_inputs$new_data, assess_inputs$new_parameters)
    #mod_out$opt$par["ln_mean_rec"]

    om_ssb <- apply(mse_tier3$naa[1:y,,1,1,s]*sable_om$dem_params$mat[1:y,,1,1]*sable_om$dem_params$waa[1:y,,1,1], 1, sum)

    assess_ssb <- SpatialSablefishAssessment::get_SSB(mod_out$report)

    assess_ssb$om_ssb <- om_ssb
    assess_ssb$sim <- s
    assess_ssb$om <- om
    return(assess_ssb)
}