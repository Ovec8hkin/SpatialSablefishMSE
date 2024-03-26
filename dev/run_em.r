rm(list=ls())
library(SpatialSablefishAssessment)

om <- readRDS("data/om2.RDS")


load(system.file("testdata", "MockAssessmentModel.RData",package="SpatialSablefishAssessment"))
nyears <- 150
extra_years <- nyears-63

# extract ages
ages <- as.numeric(colnames(om$dem_params$waa))
data$ages <- ages

# extract years
years <- 1960-1+as.numeric(rownames(om$dem_params$waa[1:nyears,,,]))
data$years <- years

# extract maturity
maturity <- t(om$dem_params$mat[1:nyears,,1,1])
data$maturity <- maturity

# extract mortality
mortality <- t(om$dem_params$mort[1:nyears,,1,1])
data$M <- mortality

# extract WAA
waa <- om$dem_params$waa[1:nyears,,,1, drop=FALSE]
waa_m <- t(waa[,,2,])
waa_f <- t(waa[,,1,])

data$male_mean_weight_by_age <- waa_m
data$female_mean_weight_by_age <- waa_f

# age-length transition matrix
male_al_trans_mat <- extend_3darray_last_dim(data$male_age_length_transition, n=extra_years)
female_al_trans_mat <- extend_3darray_last_dim(data$female_age_length_transition, n=extra_years)

data$male_age_length_transition <- male_al_trans_mat
data$female_age_length_transition <- female_al_trans_mat

# proportion male
prop_male <- rep(0.5, nyears)
data$proportion_male <- prop_male
data$proportion_male2 <- prop_male

# spawning time proportion
spawn_time <- rep(0, nyears)
data$spawning_time_proportion <- spawn_time

# fishery catch vectors
landed_caa <- subset_matrix(om$land_caa[1:nyears,,,,,1,drop=FALSE], 1, 6)
ll_catch <- apply(landed_caa[1:nyears,,,,1], 1, sum, na.rm=TRUE)
tw_catch <- apply(landed_caa[1:nyears,,,,2], 1, sum, na.rm=TRUE)

data$ll_fishery_catch <- ll_catch
data$trwl_fishery_catch <- tw_catch

# fishery selectivity
ll_sel_yr_indic <- extend_vec_last_val(data$ll_sel_by_year_indicator, n=extra_years)
tw_sel_yr_indic <- extend_vec_last_val(data$trwl_sel_by_year_indicator, n=extra_years)
jp_sel_yr_indic <- extend_vec_last_val(data$srv_jap_fishery_ll_sel_by_year_indicator, n=extra_years)

data$ll_sel_by_year_indicator <- ll_sel_yr_indic
data$trwl_sel_by_year_indicator <- tw_sel_yr_indic
data$srv_jap_fishery_ll_sel_by_year_indicator <- jp_sel_yr_indic

# survey selectivity and catchability
srv_ll_sel_yr_indic <- extend_vec_last_val(data$srv_dom_ll_sel_by_year_indicator, n=extra_years)
srv_tw_sel_yr_indic <- extend_vec_last_val(data$srv_nmfs_trwl_sel_by_year_indicator, n=extra_years)
srv_jp_sel_yr_indic <- extend_vec_last_val(data$srv_jap_ll_sel_by_year_indicator, n=extra_years)

data$srv_dom_ll_sel_by_year_indicator <- srv_ll_sel_yr_indic
data$srv_nmfs_trwl_sel_by_year_indicator <- srv_tw_sel_yr_indic
data$srv_jap_ll_sel_by_year_indicator <- srv_jp_sel_yr_indic

srv_ll_q_yr_indic <- extend_vec_last_val(data$srv_dom_ll_q_by_year_indicator, n=extra_years)
srv_tw_q_yr_indic <- extend_vec_last_val(data$srv_nmfs_trwl_q_by_year_indicator, n=extra_years)
srv_jp_q_yr_indic <- extend_vec_last_val(data$srv_jap_ll_q_by_year_indicator, n=extra_years)
jp_q_yr_indic <- extend_vec_last_val(data$srv_jap_fishery_ll_q_by_year_indicator, n=extra_years)

data$srv_dom_ll_q_by_year_indicator <- srv_ll_q_yr_indic
data$srv_nmfs_trwl_q_by_year_indicator <- srv_tw_q_yr_indic
data$srv_jap_ll_q_by_year_indicator <- srv_jp_q_yr_indic
data$srv_jap_fishery_ll_q_by_year_indicator <- jp_q_yr_indic

ll_cpue_q_yr_indic <- extend_vec_last_val(data$ll_cpue_q_by_year_indicator, n=extra_years)
data$ll_cpue_q_by_year_indicator <- ll_cpue_q_yr_indic

# catch at age data
ll_caa_indic <- extend_vec_last_val(data$ll_catchatage_indicator, n=extra_years)
ll_caa_indic[63:(length(ll_caa_indic)-1)] <- 1
data$ll_catchatage_indicator <- ll_caa_indic

obs_ll_caa <- apply(om$survey_obs$fxfish_acs[which(ll_caa_indic == 1),,,,drop=FALSE], c(1, 2), sum)
data$obs_ll_catchatage <- apply(obs_ll_caa, 1, as.double)

srv_ll_caa_indic <- extend_vec_last_val(data$srv_dom_ll_age_indicator, n=extra_years)
srv_ll_caa_indic[63:(length(years)-1)] <- 1
data$srv_dom_ll_age_indicator <- srv_ll_caa_indic

obs_srv_ll_caa <- apply(om$survey_obs$ll_acs[which(srv_ll_caa_indic == 1),,,,drop=FALSE], c(1, 2), sum)
data$obs_srv_dom_ll_age <- apply(obs_srv_ll_caa, 1, as.double)

srv_jpll_caa_indic <- extend_vec_last_val(data$srv_jap_ll_age_indicator, n=extra_years)
data$srv_jap_ll_age_indicator <- srv_jpll_caa_indic

srv_tw_caa_indic <- c(data$srv_nmfs_trwl_age_indicator, rep(0, length.out=extra_years))
data$srv_nmfs_trwl_age_indicator <- srv_tw_caa_indic

# catch at length
ll_cal_indic <- extend_vec_last_val(data$ll_catchatlgth_indicator, n=extra_years)
tw_cal_indic <- extend_vec_last_val(data$trwl_catchatlgth_indicator, n=extra_years)

data$ll_catchatlgth_indicator <- ll_cal_indic
data$trwl_catchatlgth_indicator <- tw_cal_indic

srv_ll_cal_indic <- c(data$srv_dom_ll_lgth_indicator, rep(0, length.out=extra_years))
srv_tw_cal_indic <- c(data$srv_nmfs_trwl_lgth_indicator, rep(0, length.out=extra_years))

data$srv_dom_ll_lgth_indicator <- srv_ll_cal_indic
data$srv_nmfs_trwl_lgth_indicator <- srv_tw_cal_indic

srv_jpll_cal_indic <- extend_vec_last_val(data$srv_jap_ll_lgth_indicator, n=extra_years)
jp_ll_cal_indic <- extend_vec_last_val(data$srv_jap_fishery_ll_lgth_indicator, n=extra_years)

data$srv_jap_ll_lgth_indicator <- srv_jpll_cal_indic
data$srv_jap_fishery_ll_lgth_indicator <- jp_ll_cal_indic

# Survey data
srv_ll_rpn_indic <- extend_vec_last_val(data$srv_dom_ll_bio_indicator, n=extra_years)
srv_ll_rpn_obs <- om$survey_obs$ll_rpn[srv_ll_rpn_indic]
srv_ll_rpn_ses <- om$model_options$obs_pars$surv_ll$rpn_cv*om$survey_obs$ll_rpn[srv_ll_rpn_indic]

data$srv_dom_ll_bio_indicator <- srv_ll_rpn_indic
data$obs_dom_ll_bio <- srv_ll_rpn_obs
data$se_dom_ll_bio <- srv_ll_rpn_ses

srv_jpll_indic <- extend_vec_last_val(data$srv_jap_ll_bio_indicator, n=extra_years)
data$srv_jap_ll_bio_indicator <- srv_jpll_indic

srv_tw_rpw_indic <- c(data$srv_nmfs_trwl_bio_indicator, rep(c(0, 1), length.out=extra_years))
srv_tw_rpw_obs <- om$survey_obs$tw_rpw[srv_tw_rpw_indic]
srv_tw_rpw_ses <- om$model_options$obs_pars$surv_tw$rpw_cv*om$survey_obs$tw_rpw[srv_tw_rpw_indic]

data$srv_nmfs_trwl_bio_indicator <- srv_tw_rpw_indic
data$obs_nmfs_trwl_bio <- srv_tw_rpw_obs
data$se_nmfs_trwl_bio <- srv_tw_rpw_ses

ll_cpue_indic <- extend_vec_last_val(data$ll_cpue_indicator, n=extra_years)
jp_ll_rpw_indic <- extend_vec_last_val(data$srv_jap_fishery_ll_bio_indicator, n=extra_years)

data$ll_cpue_indicator <- ll_cpue_indic
data$srv_jap_fishery_ll_bio_indicator <- jp_ll_rpw_indic

# Parameters

ln_M_year_devs <- extend_vec_last_val(parameters$ln_M_year_devs, n=extra_years)
ln_M_age_devs  <- extend_vec_last_val(parameters$ln_M_age_devs, n=extra_years)
ln_rec_dev <- c(parameters$ln_rec_dev, rep(0, extra_years))
ln_ll_F <- c(parameters$ln_ll_F_devs, rep(0, extra_years))
ln_tw_F <- c(parameters$ln_trwl_F_devs, rep(0, extra_years))

parameters$ln_M_year_devs <- ln_M_year_devs
parameters$ln_M_age_devs <- ln_M_age_devs
parameters$ln_rec_dev <- ln_rec_dev
parameters$ln_ll_F_devs <- ln_ll_F
parameters$ln_trwl_F_devs <- ln_tw_F

validate_input_data_and_parameters(data, parameters)

my_model = TMB::MakeADFun(data = data,
                               parameters = parameters,
                               DLL = "SpatialSablefishAssessment_TMBExports")

mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))

mle_report = my_model$report(mle_optim$par)  

true_ssb <- apply(om$naa[1:nyears,,1,1,1]*om$dem_params$mat[1:nyears,,1,1]*om$dem_params$waa[1:nyears,,1,1], 1, sum)

ssb <- get_SSB(mle_report) %>% mutate(om=true_ssb)
ggplot(ssb)+
    geom_line(aes(x=Year, y=SSB))+
    geom_line(aes(x=Year, y=om), color="red")+
    scale_y_continuous(limits=c(0, 300))+
    theme_bw()
