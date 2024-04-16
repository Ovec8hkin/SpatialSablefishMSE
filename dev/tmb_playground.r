#' Playground space to investigate problems with the TMB Estimation model
#' 
#' This script can be used to investigate alternative parameterizations
#' of the TMB estimation model being used as the EM for the Alaska
#' sablefish MSE. The TMB model estimated of SSB and LL fishing mortality
#' are known to be biased due to issues with estimating selectivity.
#' 
#' @author Joshua A. Zahner ('jzahner@alaska.edu')
#' @author Matt L. Cheng
#' 

library(TMB)
library(here)
library(tidyverse)
library(ggdist)
library(doParallel)
library(patchwork)
library(cowplot)

#setwd("~/Desktop/Side projects/SablefishMSE")
source("R/run_mse_multiple.R")
source("R/run_mse.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/data_utils.R")
source("R/format_em_data.R")
source("R/fit_TMB_model.R")
source("R/reference_points.R")

library(devtools)
devtools::load_all("~/Desktop/Projects/afscOM")


sable_om <- readRDS("data/sablefish_om.RDS") # Read this saved OM from a file
sable_om$model_options$simulate_observations <- TRUE # Enable simulating observations
sable_om$model_options$obs_pars$surv_ll$as_integers = TRUE
sable_om$model_options$obs_pars$surv_tw$as_integers = TRUE
sable_om$model_options$obs_pars$fish_fx$as_integers = TRUE
sable_om$model_options$run_estimation = FALSE

tier3 <- function(ref_pts, naa, dem_params){
  ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
  return(
    npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
  )
}

nsims <- 100
set.seed(1120)
seeds <- sample(1:1e3, nsims)
nyears <- 100
sable_om$model_options$obs_pars$fish_fx$ac_samps <- 30
sable_om$model_options$obs_pars$surv_ll$ac_samps <- 30
sable_om$model_options$obs_pars$surv_tw$ac_samps <- 30
sable_om$model_options$obs_pars$surv_ll$rpn_cv <- 0.1
sable_om$model_options$obs_pars$surv_tw$rpw_cv <- 0.1

s <- 1
mse_tier3 <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=nyears, om=sable_om, hcr=tier3)

data <- data.frame()

cores <- parallel::detectCores()-2
cl <- parallel::makeCluster(cores, outfile="")
registerDoParallel(cl)

out <- pbapply::pblapply(1:nsims, function(s, sable_om, mse_tier3, nyears){
    library(tidyverse)
    library(TMB)

    source("R/format_em_data.R")
    source("R/fit_TMB_model.R")
    source("R/reference_points.R")
    
    library(devtools)
    devtools::load_all("~/Desktop/Projects/afscOM")

    i=nyears-63
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
    
    #' Fit the model
    mod_out <- fit_TMB_model(assess_inputs$new_data, assess_inputs$new_parameters)
    
    mle_report = mod_out$report
    mle_optim = mod_out$opt

    est_rec <- SpatialSablefishAssessment::get_recruitment(mle_report)
    om_rec <- apply(mse_tier3$naa[1:nyears,1,,,s], 1, sum)
    
    est_ssb <- SpatialSablefishAssessment::get_SSB(mle_report)
    om_ssb <- apply(mse_tier3$naa[1:nyears,,1,1,s]*sable_om$dem_params$mat[1:nyears,,1,1]*sable_om$dem_params$waa[1:nyears,,1,1], 1, sum)
    
    est_Fs <- SpatialSablefishAssessment::get_fishing_mortalities(mle_report) %>% as_tibble() %>% group_by(Year) %>% summarise(F = sum(F))
    om_Fs <- mse_tier3$out_f[1:nyears,1, 1, 1, s]
    
    est_fish_Fs <- SpatialSablefishAssessment::get_fishing_mortalities(mle_report) %>% as_tibble() %>%
        mutate(
            true_Fs = c(apply(mse_tier3$faa[1:nyears,,,1,1,s], 1, max), apply(mse_tier3$faa[1:nyears,,,1,2,s], 1, max)),
            avg_F = c(rep(exp(mle_optim$par["ln_ll_F_avg"]), nyears), rep(exp(mle_optim$par["ln_trwl_F_avg"]), nyears)),
            true_avg_F = c(rep(exp(NA), nyears), rep(exp(NA), nyears))
            )
    
    est_catches <- SpatialSablefishAssessment::get_catches(mle_report) %>% as_tibble() %>%
        pivot_wider(names_from=type, values_from=Catch)
    
    selex <- SpatialSablefishAssessment::get_selectivities(mle_report)

    sel_est <- array(NA, dim=c(1, 30, 1, 1, 2))
    sel_est[1,,1,1,1] <- selex %>% filter(gear == "fixed", sex == "female", time_block == 3) %>% pull(value)
    sel_est[1,,1,1,2] <- selex %>% filter(gear == "trawl", sex == "female", time_block == 1) %>% pull(value)

    joint_self_true <- apply(sable_om$dem_params$sel[nyears,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(sable_om$dem_params$sel[nyears,,1,,,drop=FALSE], c(1, 2), sum))
    joint_self_est <- apply(sel_est, c(1, 2), sum)/max(apply(sel_est, c(1, 2), sum))
    rps_est <-  calculate_ref_points(30, mle_report$M[1,1], sable_om$dem_params$mat[nyears,,1,1], sable_om$dem_params$waa[nyears,,1,1], joint_self_est, rep(1, 30), avg_rec=mean(est_rec$Recruitment[17:nyears]/2))
    rps_true <- calculate_ref_points(30, sable_om$dem_params$mort[1,1,1,1], sable_om$dem_params$mat[nyears,,1,1], sable_om$dem_params$waa[nyears,,1,1], joint_self_true, rep(1, 30), avg_rec=mean(mse_tier3$naa[17:nyears, 1, 1, 1, s]))

    tmp <- data.frame(
        Year=1960:(1960+nyears-1), 
        est_rec=est_rec$Recruitment, 
        om_rec=om_rec, 
        est_ssb=est_ssb$SSB, 
        om_ssb=om_ssb, 
        est_F=est_Fs$F, 
        om_F=om_Fs, 
        est_ll_f=est_fish_Fs %>% filter(Fishery == "Fixed gear") %>% pull(F),
        om_ll_f=est_fish_Fs %>% filter(Fishery == "Fixed gear") %>% pull(true_Fs),
        est_tw_f=est_fish_Fs %>% filter(Fishery == "Trawl gear") %>% pull(F),
        om_tw_f=est_fish_Fs %>% filter(Fishery == "Trawl gear") %>% pull(true_Fs),
        est_ll_catch = est_catches %>% filter(Fishery == "Fixed gear") %>% pull(Predicted),
        om_ll_catch = est_catches %>% filter(Fishery == "Fixed gear") %>% pull(Observed),
        est_tw_catch = est_catches %>% filter(Fishery == "Trawl") %>% pull(Predicted),
        om_tw_catch = est_catches %>% filter(Fishery == "Trawl") %>% pull(Observed),
        # B0 = mle_report$Bzero,
        # R0 = mle_report$mean_rec/2,
        F40_est = rps_est$F40,
        F35_est = rps_est$F35,
        B40_est = rps_est$B40,
        F40_true = rps_true$F40,
        F35_true = rps_true$F35,
        B40_true = rps_true$B40,
        sim=s
    )

    return(tmp)
}, sable_om=sable_om, mse_tier3=mse_tier3, nyears=nyears, cl=cl)

data <- bind_rows(out)

stopCluster(cl)

df <- data %>% group_by(Year) %>%
  ggdist::median_qi(est_rec, om_rec, est_ssb, om_ssb, est_F, 
                    om_F, est_ll_f, om_ll_f, est_tw_f, om_tw_f, 
                    est_ll_catch, om_ll_catch, est_tw_catch, om_tw_catch, .width=c(0.50, 0.80))

rec_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_rec, ymin=om_rec.lower, ymax=om_rec.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_rec, ymin=est_rec.lower, ymax=est_rec.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 200))+
  ggtitle("Recruitment")+
  guides(fill="none")+
  theme_bw()

ssb_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_ssb, ymin=om_ssb.lower, ymax=om_ssb.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_ssb, ymin=est_ssb.lower, ymax=est_ssb.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 350))+
  ggtitle("Spawning Biomass")+
  guides(fill="none")+
  theme_bw()

F_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_F, ymin=om_F.lower, ymax=om_F.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_F, ymin=est_F.lower, ymax=est_F.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 0.3))+
  ggtitle("Fishing Mortality")+
  guides(fill="none")+
  theme_bw()

ll_F_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_ll_f, ymin=om_ll_f.lower, ymax=om_ll_f.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_ll_f, ymin=est_ll_f.lower, ymax=est_ll_f.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 0.15))+
  ggtitle("Fixed Gear Fishing Mortality")+
  guides(fill="none")+
  theme_bw()

tw_F_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_tw_f, ymin=om_tw_f.lower, ymax=om_tw_f.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_tw_f, ymin=est_tw_f.lower, ymax=est_tw_f.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 0.12))+
  ggtitle("Trawl Gear Fishing Mortality")+
  guides(fill="none")+
  theme_bw()

ll_catch_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_ll_catch, ymin=om_ll_catch.lower, ymax=om_ll_catch.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_ll_catch, ymin=est_ll_catch.lower, ymax=est_ll_catch.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 50))+
  ggtitle("Fixed Gear Catch")+
  guides(fill="none")+
  theme_bw()

tw_catch_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_tw_catch, ymin=om_tw_catch.lower, ymax=om_tw_catch.upper), size=0.6)+
  geom_line(aes(x=Year, y=est_tw_catch, ymin=est_tw_catch.lower, ymax=est_tw_catch.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 50))+
  ggtitle("Trawl Gear Catch")+
  guides(fill="none")+
  theme_bw()

plot_grid(
  plot_grid(rec_plot, ssb_plot, F_plot, ncol=3),
  plot_grid(
    plot_grid(ll_F_plot, tw_F_plot, ncol=1), 
    plot_grid(ll_catch_plot, tw_catch_plot, ncol=1), 
    ncol=2
  ),
  ncol=1
)

#' Plot reference points (F35 (F_OFL), F40 (F_ABC), and B40) to look for
#' estimation bias.
data %>% select(F40_est, F40_true, F35_est, F35_true, B40_est, B40_true, sim) %>%
  median_qi(sim, F40_est, F40_true, F35_est, F35_true, B40_est, B40_true, .width=c(0.50, 0.95)) %>%
  reformat_ggdist_long(1) %>%
  separate(name, into=c("RP", "model"), sep="_") %>%

  ggplot()+
    geom_pointrange(aes(x=model, y=median, ymin=lower, ymax=upper, color="model"))+
    facet_wrap(~RP, scales="free_y")+
    ggh4x::facetted_pos_scales(
      y = list(
        scale_y_continuous(limits=c(50, 150)),
        scale_y_continuous(limits=c(0.1, 0.15)),
        scale_y_continuous(limits=c(0.09, 0.11))
      )
    )+
    theme_bw()


#' Plot selectivity curves to evaluate fits to OM selectivity
s <- 1
i <- nyears-63
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

#' Fit the model
mod_out <- fit_TMB_model(assess_inputs$new_data, assess_inputs$new_parameters)
mle_report = mod_out$report


selex <- SpatialSablefishAssessment::get_selectivities(mle_report)
om_selex <- reshape2::melt(sable_om$dem_params$sel) %>% as_tibble %>%
    filter(time %in% c(1, 45, 65)) %>%
    select(-c(region)) %>%
    mutate(
        time = case_when(time == 1 ~ 1, time == 45 ~ 2, time == 65 ~ 3),
        sex = case_when(sex == "F" ~ "female", sex == "M" ~ "male"),
        fleet = case_when(fleet == "Fixed" ~ "fixed", fleet == "Trawl" ~ "trawl"),
        age = age-1
    ) %>%
    bind_rows(
        reshape2::melt(sable_om$dem_params$surv_sel) %>% as_tibble %>%
            filter(time %in% c(1, 65)) %>%
            select(-c(region)) %>%
            mutate(
                time = case_when(time == 1 ~ 1, time == 65 ~ 2),
                sex = case_when(sex == "F" ~ "female", sex == "M" ~ "male"),
                fleet = case_when(fleet == "Fixed" ~ "domsurveyll", fleet == "Trawl" ~ "nmfssurveytrwl"),
                age = age-1
            )

    ) %>%
    rename("time_block"=time, "gear"=fleet) %>%
    arrange(age, value, gear, sex, time_block) %>%
    print(n=10)

selex %>% left_join(om_selex, by=c("time_block", "age", "sex", "gear")) %>%
    rename("est"=value.x, "true"=value.y) %>%
    filter(sex != "combined", gear != "japsurveyll") %>%
    ggplot()+
        geom_line(aes(x=age, y=est, color=time_block, group=time_block))+
        geom_line(aes(x=age, y=true, color=time_block, group=time_block), linetype="dashed")+
        facet_grid(vars(gear), vars(sex))+
        theme_bw()




a
