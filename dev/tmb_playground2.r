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
rm(list=ls())

library(TMB)
library(here)
library(tidyverse)
library(ggdist)
library(doParallel)
library(patchwork)
library(cowplot)
library(ggh4x)

#setwd("~/Desktop/Side projects/SablefishMSE")
source("R/run_mse_parallel.R")
source("R/run_mse.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/data_utils.R")
source("R/format_em_data.R")
source("R/fit_TMB_model.R")
source("R/reference_points.R")
source("R/recruitment_utils.R")

library(devtools)
devtools::load_all("~/Desktop/Projects/afscOM")

nyears <- 165

sable_om <- readRDS("data/sablefish_om_big.RDS") # Read this saved OM from a file

tier3 <- function(ref_pts, naa, dem_params){
  ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
  return(
    npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
  )
}

assessment <- dget("data/sablefish_assessment_2023.rdat")
hist_recruits <- assessment$natage.female[,1]*2

sable_om$model_options$recruitment$func <- resample_recruits
sable_om$model_options$recruitment$pars <- list(
    hist_recruits = hist_recruits,
    nyears = nyears - length(hist_recruits) + 1
)

nsims <- 20
set.seed(94265)
seeds <- sample(1:1e3, nsims)

s <- 1
mse_tier3 <- run_mse_parallel(nsims=nsims, seeds=seeds, nyears=nyears, om=sable_om, hcr=tier3)

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

    i <- nyears-63
    y <- 63+i
    assess_inputs <- simulate_em_data_sex_disaggregate(
        nyears = y,
        dem_params = afscOM::subset_dem_params(sable_om$dem_params, 1:y, d=1, drop=FALSE),
        land_caa = afscOM::subset_matrix(mse_tier3$land_caa[1:y,,,,,s,drop=FALSE], 1, d=6, drop=TRUE),
        survey_indices = afscOM::subset_dem_params(afscOM::subset_dem_params(mse_tier3$survey_obs, 1:y, d=1, drop=FALSE), s, d=5, drop=TRUE),
        fxfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$fxfish_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
        twfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$twfish_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
        ll_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$ll_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
        tw_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$tw_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
        model_options = sable_om$model_options,
        added_years = i,
        file_suffix = s
    )

    # assess_inputs$new_data$mu_srv_dom_ll_q <- 6.412379
    # assess_inputs$new_data$sd_srv_dom_ll_q <- 0.05
    # assess_inputs$new_data$mu_srv_nmfs_trwl_q <- 0.8580096
    # assess_inputs$new_data$sd_srv_nmfs_trwl_q <- 0.05
    # assess_inputs$new_data$loglik_wgt_q_priors <- 1

    # assess_inputs$new_data$mu_M <- sable_om$dem_params$mort[1,1,1,1]
    # assess_inputs$new_data$sd_M <- 8.5481e-002

    # fix_pars = list(
    #   ln_M = log(sable_om$dem_params$mort[1,1,1,1])
    #   # ln_srv_nmfs_trwl_q = -0.153144044674,
    #   # ln_srv_dom_ll_q = 1.85823035364
    # )

    fix_pars <- NULL

    mod_out <- fit_TMB_model(
        data = assess_inputs$new_data,
        parameters = assess_inputs$new_parameters,
        model_name = "CurrentAssessmentDisaggregated",
        fix_pars = fix_pars
    )


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
    
    selex <- SpatialSablefishAssessment::get_selectivities(mle_report) %>% mutate(sim=s)

    sel_est <- array(NA, dim=c(1, 30, 1, 1, 2))
    sel_est[1,,1,1,1] <- selex %>% filter(gear == "fixed", sex == "female", time_block == 3) %>% pull(value)
    sel_est[1,,1,1,2] <- selex %>% filter(gear == "trawl", sex == "female", time_block == 1) %>% pull(value)

    prop_fs <- apply(mse_tier3$faa[y,,,1,,s, drop=FALSE], 5, max)/sum(apply(mse_tier3$faa[y,,,1,,s, drop=FALSE], 5, max))
    joint_self_true <- apply(sable_om$dem_params$sel[nyears,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum)/max(apply(sable_om$dem_params$sel[nyears,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum))
    
    spinup_years <- 63
    F_ll_f <- t(mle_report$F_ll_f[,y])
    F_ll_m <- t(mle_report$F_ll_m[,y])
    F_tw_f <- t(mle_report$F_trwl_f[,y])
    F_tw_m <- t(mle_report$F_trwl_m[,y])
    faa_est <- array(NA, dim=c(1, 30, 2, 1, 2))
    faa_est[1,,1,,1] <- F_ll_f
    faa_est[1,,2,,1] <- F_ll_m
    faa_est[1,,1,,2] <- F_tw_f
    faa_est[1,,2,,2] <- F_tw_m
    prop_fs <- apply(faa_est[1,,,1,, drop=FALSE], 5, max)/sum(apply(faa_est[1,,,1,, drop=FALSE], 5, max))
    joint_self_est <- apply(sel_est*prop_fs, c(1, 2), sum)/max(apply(sel_est*prop_fs, c(1, 2), sum))
    
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
        est_F40 = rps_est$F40,
        est_F35 = rps_est$F35,
        est_B40 = rps_est$B40,
        om_F40 = rps_true$F40,
        om_F35 = rps_true$F35,
        om_B40 = rps_true$B40,
        M = exp(mle_optim$par["ln_M"]),
        ll_q = exp(mle_optim$par["ln_srv_dom_ll_q"]),
        tw_q = exp(mle_optim$par["ln_srv_nmfs_trwl_q"]),
        sim=s
    )

    return(afscOM::listN(tmp, selex))
}, sable_om=sable_om, mse_tier3=mse_tier3, nyears=nyears, cl=cl)

data <- bind_rows(lapply(out, \(x) x$tmp))

stopCluster(cl)

summary(data[, c("ll_q", "tw_q", "M")])

df <- data %>% group_by(Year) %>%
  ggdist::median_qi(est_rec, om_rec, est_ssb, om_ssb, est_F, 
                    om_F, est_ll_f, om_ll_f, est_tw_f, om_tw_f, 
                    est_ll_catch, om_ll_catch, est_tw_catch, om_tw_catch, .width=c(0.50, 0.80))

rec_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_rec, ymin=om_rec.lower, ymax=om_rec.upper), size=0.6)+
#   geom_pointrange(aes(x=Year, y=est_rec, ymin=est_rec.lower, ymax=est_rec.upper), color="red", size=0.1)+
  geom_line(aes(x=Year, y=est_rec), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 200))+
  ggtitle("Recruitment")+
  guides(fill="none")+
  theme_bw()

ssb_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_ssb, ymin=om_ssb.lower, ymax=om_ssb.upper), size=0.6)+
#   geom_pointrange(aes(x=Year, y=est_ssb, ymin=est_ssb.lower, ymax=est_ssb.upper), color="red", size=0.1)+
  geom_line(aes(x=Year, y=est_ssb), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 350))+
  ggtitle("Spawning Biomass")+
  guides(fill="none")+
  theme_bw()

F_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_F, ymin=om_F.lower, ymax=om_F.upper), size=0.6)+
#   geom_pointrange(aes(x=Year, y=est_F, ymin=est_F.lower, ymax=est_F.upper), color="red", size=0.1)+
  geom_line(aes(x=Year, y=est_F, ymin=est_F.lower, ymax=est_F.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 0.3))+
  ggtitle("Fishing Mortality")+
  guides(fill="none")+
  theme_bw()

ll_F_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_ll_f, ymin=om_ll_f.lower, ymax=om_ll_f.upper), size=0.6)+
#   geom_pointrange(aes(x=Year, y=est_ll_f, ymin=est_ll_f.lower, ymax=est_ll_f.upper), color="red", size=0.1)+
  geom_line(aes(x=Year, y=est_ll_f, ymin=est_ll_f.lower, ymax=est_ll_f.upper), color="red")+
  scale_fill_brewer()+
  scale_y_continuous(limits=c(0, 0.15))+
  ggtitle("Fixed Gear Fishing Mortality")+
  guides(fill="none")+
  theme_bw()

tw_F_plot <- ggplot(df)+
  geom_lineribbon(aes(x=Year, y=om_tw_f, ymin=om_tw_f.lower, ymax=om_tw_f.upper), size=0.6)+
#   geom_pointrange(aes(x=Year, y=est_tw_f, ymin=est_tw_f.lower, ymax=est_tw_f.upper), color="red", size=0.1)+
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

summary_plot <- plot_grid(
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
ref_pts <- data %>% select(contains("40"), contains("35"), sim) %>%
  as_tibble() %>%
  median_qi(sim, est_F40, om_F40, est_F35, om_F35, est_B40, om_B40, .width=c(0.50, 0.95)) %>%
  reformat_ggdist_long(1) %>%
  separate(name, into=c("model", "RP"), sep="_")

ref_pts_plot <- ggplot(ref_pts)+
  geom_pointrange(aes(x=model, y=median, ymin=lower, ymax=upper, color=model))+
  facet_wrap(~RP, scales="free_y")+
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits=c(75, 125)),
      scale_y_continuous(limits=c(0.15, 0.2)),
      scale_y_continuous(limits=c(0.13, 0.18))
    )
  )+
  theme_bw()


#' Plot selectivity curves to evaluate fits to OM selectivity
selex <- bind_rows(lapply(out, \(x) x$selex))
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

selectivities <- selex %>% left_join(om_selex, by=c("time_block", "age", "sex", "gear")) %>%
    rename("est"=value.x, "true"=value.y) %>%
    filter(sex != "combined", gear != "japsurveyll") %>%
    group_by(age, gear, sex, time_block) %>%
    median_qi(est, true, .width=c(0.50, 0.80)) %>%
    reformat_ggdist_long(n=4) %>%
    mutate(
      time_block = factor(time_block)
    )

selex_plot <- ggplot(selectivities)+
    geom_lineribbon(aes(x=age, y=median, ymin=lower, ymax=upper, color=time_block, group=time_block, linetype=name), size=0.6)+
    scale_color_manual(values=c("black", "red", "blue")) + 
    scale_fill_brewer()+
    facet_grid(vars(gear), vars(sex))+
    theme_bw()

#' Make relative error plots
rel_error <- data %>% as_tibble() %>%
  pivot_longer(-c(Year, sim), values_to="value", names_to="quantity") %>%
  filter(quantity != "M") %>%
  separate(quantity, c("model", "quantity"), extra = "merge") %>%
  pivot_wider(names_from="model", values_from="value") %>%
  mutate(
    rel_error = (est-om)/om
  ) %>%
  select(Year, sim, quantity, rel_error) %>%
  pivot_wider(names_from="quantity", values_from="rel_error") %>%
  group_by(Year) %>%
  median_qi(rec, ssb, F, ll_f, tw_f, B40, F40, F35, .width=c(0.50, 0.80)) %>%
  reformat_ggdist_long(1) %>%
  mutate(
    name = factor(name, levels=c("ssb", "F", "rec", "ll_f", "tw_f", "F35", "F40", "B40"))
  )

avg_rel_error <- rel_error %>% 
  group_by(name) %>% 
  summarise(avg = median(median[Year > 2022])) %>%
  mutate(
    text_y = case_when(
      name == "ssb" ~ 0.085,
      name == "F" ~ 0.30,
      name == "rec" ~ 0.375,
      name == "ll_f" ~ 0.30,
      name == "tw_f" ~ 0.30,
      name == "F35" ~ 0.04,
      name == "F40" ~ 0.04,
      name == "B40" ~ 0.04,
    )
  )

rel_error_plot <- ggplot(rel_error)+
  geom_lineribbon(aes(x=Year, y=median, ymin=lower, ymax=upper), size=0.6)+
  geom_vline(xintercept=2022)+
  geom_hline(data=avg_rel_error, aes(yintercept=avg), linetype="solid", color="red", size=0.6) + 
  geom_hline(yintercept=0.0, size=1)+
  geom_text(data=avg_rel_error, aes(x=2100, y=text_y, label=paste0(round(100*avg, 3), "%")), color="red")+
  scale_fill_brewer()+
  facet_wrap(~name, scales="free_y")+
  facetted_pos_scales(
      y=list(
          scale_y_continuous(limits=c(-0.12, 0.12)),
          scale_y_continuous(limits=c(-0.05, 0.35)),
          scale_y_continuous(limits=c(-0.5, 0.5)),
          scale_y_continuous(limits=c(-0.05, 0.35)),
          scale_y_continuous(limits=c(-0.05, 0.35)),
          scale_y_continuous(limits=c(-0.05, 0.05)),
          scale_y_continuous(limits=c(-0.05, 0.05)),
          scale_y_continuous(limits=c(-0.05, 0.05))
      )
  )+
  theme_bw()


summary_plot
rel_error_plot
# ref_pts_plot
# selex_plot


a

########

s <- 1
i <- nyears-63
y <- 63+i
assess_inputs <- simulate_em_data_sex_disaggregate(
    nyears = y,
    dem_params = afscOM::subset_dem_params(sable_om$dem_params, 1:y, d=1, drop=FALSE),
    land_caa = afscOM::subset_matrix(mse_tier3$land_caa[1:y,,,,,s,drop=FALSE], 1, d=6, drop=TRUE),
    survey_indices = afscOM::subset_dem_params(afscOM::subset_dem_params(mse_tier3$survey_obs, 1:y, d=1, drop=FALSE), s, d=5, drop=TRUE),
    fxfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$fxfish_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
    twfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$twfish_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
    ll_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$ll_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
    tw_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$tw_acs[1:y,,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
    model_options = sable_om$model_options,
    added_years = i,
    file_suffix = s
)

vars <- c(
  "ll_catch",
  "trwl_catch",
  "ll_catchatage_indicator",
  "ll_catchatlgth_indicator",
  "trwl_catchatlgth_indicator",
  "srv_dom_ll_bio_indicator",
  "srv_jap_ll_bio_indicator",
  "ll_cpue_indicator",
  "srv_dom_ll_age_indicator",
  "srv_dom_ll_lgth_indicator",
  "srv_jap_ll_age_indicator",
  "srv_jap_ll_lgth_indicator",
  "srv_nmfs_trwl_age_indicator",
  "srv_nmfs_trwl_lgth_indicator",
  "srv_nmfs_trwl_bio_indicator",
  "srv_jap_fishery_ll_bio_indicator",
  "srv_jap_fishery_ll_lgth_indicator",
  "trwl_catchatage_indicator"
)

data_df <- data.frame(Year=1960:(1960+nyears-1))
for(v in vars){
  if(v == "ll_catch" || v == "trwl_catch"){
    data_df[,v] <- 1
  }else{
    data_df[,v] <- assess_inputs$new_data[v]
  }
}

data_df %>% as_tibble() %>%
  pivot_longer(-c("Year"), names_to="dataset", values_to="value") %>%
  mutate(
    datatype = case_when(
        grepl("age", dataset) ~ "Age Composition",
        grepl("lgth", dataset) ~ "Length Composition",
        grepl("catch", dataset) ~ "Catches",
        grepl("bio", dataset) ~ "Survey Index",
        grepl("cpue", dataset) ~ "CPUE Index"
    ),

    fleet = case_when(
      grepl("srv_nmfs_trwl", dataset) ~ "Trawl Survey",
      grepl("srv_dom_ll", dataset) ~ "Longline Survey",
      grepl("srv_jap", dataset) ~ "Japanese LL Fishery/Survey",
      grepl("ll", dataset) ~ "Longline Fishery",
      grepl("trwl", dataset) ~ "Trawl Fishery"
    ),

    value = ifelse(grepl("srv_jap", dataset), 0, value),
    value = ifelse(grepl("cpue", dataset), 0, value),

    value = na_if(value, 0)

  ) %>%
  drop_na() %>%
  
  ggplot(aes(x=Year))+
    geom_point(aes(y=fleet, size=value/5, color=fleet))+
    scale_size_continuous(limits=c(0, 1))+
    ggforce::facet_col(~datatype, scales="free_y", space="free")+
    scale_x_continuous(limits=c(1960, 1960+nyears-1), breaks=seq(1960, 1960+nyears, 20))+
    theme_bw()+
    guides(color=guide_none(), size=guide_none())+
    theme(
      strip.background = element_blank(),
      strip.placement = "inside",
      strip.text = element_text(size=14),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(size=12),
      axis.title.y = element_blank()
    )
