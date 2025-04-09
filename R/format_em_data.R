#' Format OM Observation Data for TMB Assessment Model
#' 
#' Format data and observations from an `afscOM` operating model
#' and update the model data and parameters required by the 
#' `SpatialSablefishAssessment` TMB assessment model. 
#' 
#' This function requires the existence of 'data/sablefish_em_data_2022.RDS',
#' and 'data/sablefish_em_par_2022.RDS' files, which come packaged with the
#' `SablefishMSE` codebase. 
#'
#' @param nyears the number of years of data that will be passed to the model
#' @param dem_params demographic parameter matrices subsetted to 1 year
#' @param land_caa one year worth of landed catch-at-age data (dimensins [1, nages, nsexes, nregions, nfleets])
#' @param survey_indices one year worth of survey indices (LL RPN, LL RPW, and TW RPW)
#' @param fxfish_caa_obs one year worth of catch-at-age observation from the fixed gear fishery
#' @param ll_ac_obs one year worth of age composition observations frmo the longline survey
#' @param model_options list of model options provided to the OM
#' @param added_year the number of new years of data being added (should usually be 1) 
#'
#' @export format_em_data
#'
#' @example
#'
format_em_data <- function(nyears, dem_params, land_caa, survey_indices, fxfish_caa_obs, ll_ac_obs, model_options, added_years=1, file_suffix=""){

    # Read the most current set of data from an RDS file
    # If those files don't exist, then default back to 
    # using data through 2022.
    if(!file.exists(paste0("data/sablefish_em_data_curr_",file_suffix,".RDS"))){
        data_file <- "data/sablefish_em_data_2022.RDS"
        param_file <- "data/sablefish_em_par_2022.RDS"
    }else{
        data_file <- paste0("data/sablefish_em_data_curr_",file_suffix,".RDS")
        param_file <- paste0("data/sablefish_em_par_curr_",file_suffix,".RDS")
    }

    new_data <- readRDS(data_file)

    extra_years <- added_years

    # extract ages
    ages <- as.numeric(colnames(dem_params$waa))
    new_data$ages <- ages

    # extract years
    years <- as.double(1960:(1960+nyears-1))
    new_data$years <- years

    # extract maturity
    maturity <- dem_params$mat[,,1,1]
    if(is.null(dim(maturity))){
        maturity <- matrix(maturity, ncol=1)
    }else{
        maturity <- t(maturity)
    }
    new_data$maturity <- cbind(new_data$maturity, maturity)

    # extract mortality
    mortality <- dem_params$mort[,,1,1]
    if(is.null(dim(mortality))){
        mortality <- matrix(mortality, ncol=1)
    }else{
        mortality <- t(mortality)
    }
    new_data$M <- cbind(new_data$M, mortality)

    # extract WAA
    waa <- dem_params$waa[,,,1, drop=FALSE]
    waa_m <- waa[,,2,]
    waa_f <- waa[,,1,]

    if(is.null(dim(waa_m))){
        waa_m <- matrix(waa_m, ncol=1)
        waa_f <- matrix(waa_f, ncol=1)
    }else{
        waa_m <- t(waa_m)
        waa_f <- t(waa_f)
    }

    new_data$male_mean_weight_by_age <- cbind(new_data$male_mean_weight_by_age, waa_m)
    new_data$female_mean_weight_by_age <- cbind(new_data$female_mean_weight_by_age, waa_f)

    # age-length transition matrix
    male_al_trans_mat <- SpatialSablefishAssessment::extend_3darray_last_dim(new_data$male_age_length_transition, n=extra_years)
    female_al_trans_mat <- SpatialSablefishAssessment::extend_3darray_last_dim(new_data$female_age_length_transition, n=extra_years)

    new_data$male_age_length_transition <- male_al_trans_mat
    new_data$female_age_length_transition <- female_al_trans_mat

    # proportion male
    prop_male <- rep(0.5, nyears)
    new_data$proportion_male <- prop_male
    new_data$proportion_male2 <- prop_male

    # spawning time proportion
    spawn_time <- rep(0, nyears)
    new_data$spawning_time_proportion <- spawn_time

    # recruitment bias ramp
    new_data$end_yr_rec_est_idx <- nyears-1

    # fishery catch vectors
    #landed_caa <- subset_matrix(land_caa[,,,,,drop=FALSE], 1, 6)
    ll_catch <- apply(land_caa[,,,,1, drop=FALSE], 1, sum, na.rm=TRUE)
    tw_catch <- apply(land_caa[,,,,2, drop=FALSE], 1, sum, na.rm=TRUE)

    new_data$ll_fishery_catch <- c(new_data$ll_fishery_catch, ll_catch)
    new_data$trwl_fishery_catch <- c(new_data$trwl_fishery_catch, tw_catch)

    # fishery selectivity
    ll_sel_yr_indic     <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_sel_by_year_indicator, n=extra_years)
    tw_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$trwl_sel_by_year_indicator, n=extra_years)
    jp_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_sel_by_year_indicator, n=extra_years)

    new_data$ll_sel_by_year_indicator <- ll_sel_yr_indic
    new_data$trwl_sel_by_year_indicator <- tw_sel_yr_indic
    new_data$srv_jap_fishery_ll_sel_by_year_indicator <- jp_sel_yr_indic

    # survey selectivity and catchability
    srv_ll_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_sel_by_year_indicator, n=extra_years)
    srv_tw_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_nmfs_trwl_sel_by_year_indicator, n=extra_years)
    srv_jp_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_sel_by_year_indicator, n=extra_years)

    new_data$srv_dom_ll_sel_by_year_indicator <- srv_ll_sel_yr_indic
    new_data$srv_nmfs_trwl_sel_by_year_indicator <- srv_tw_sel_yr_indic
    new_data$srv_jap_ll_sel_by_year_indicator <- srv_jp_sel_yr_indic

    srv_ll_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_q_by_year_indicator, n=extra_years)
    srv_tw_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_nmfs_trwl_q_by_year_indicator, n=extra_years)
    srv_jp_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_q_by_year_indicator, n=extra_years)
    jp_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_q_by_year_indicator, n=extra_years)

    new_data$srv_dom_ll_q_by_year_indicator <- srv_ll_q_yr_indic
    new_data$srv_nmfs_trwl_q_by_year_indicator <- srv_tw_q_yr_indic
    new_data$srv_jap_ll_q_by_year_indicator <- srv_jp_q_yr_indic
    new_data$srv_jap_fishery_ll_q_by_year_indicator <- jp_q_yr_indic

    ll_cpue_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_cpue_q_by_year_indicator, n=extra_years)
    new_data$ll_cpue_q_by_year_indicator <- ll_cpue_q_yr_indic

    # catch at age new_data
    ll_caa_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_catchatage_indicator, n=extra_years)
    ll_caa_indic[length(new_data$ll_catchatage_indicator):(nyears-1)] <- 1
    new_data$ll_catchatage_indicator <- ll_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_ll_caa <- apply(fxfish_caa_obs[,,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_ll_catchatage <- cbind(new_data$obs_ll_catchatage, apply(obs_ll_caa, 1, as.double))

    srv_ll_caa_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_age_indicator, n=extra_years)
    srv_ll_caa_indic[length(new_data$srv_dom_ll_age_indicator):(nyears-1)] <- 1
    new_data$srv_dom_ll_age_indicator <- srv_ll_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_srv_ll_caa <- apply(ll_ac_obs[,,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_srv_dom_ll_age <- cbind(new_data$obs_srv_dom_ll_age, apply(obs_srv_ll_caa, 1, as.double))

    srv_jpll_caa_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_age_indicator, n=extra_years)
    new_data$srv_jap_ll_age_indicator <- srv_jpll_caa_indic

    srv_tw_caa_indic <- c(new_data$srv_nmfs_trwl_age_indicator, rep(0, length.out=extra_years))
    new_data$srv_nmfs_trwl_age_indicator <- srv_tw_caa_indic

    # catch at length
    ll_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_catchatlgth_indicator, n=extra_years)
    tw_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$trwl_catchatlgth_indicator, n=extra_years)

    new_data$ll_catchatlgth_indicator <- ll_cal_indic
    new_data$trwl_catchatlgth_indicator <- tw_cal_indic

    srv_ll_cal_indic <- c(new_data$srv_dom_ll_lgth_indicator, rep(0, length.out=extra_years))
    srv_tw_cal_indic <- c(new_data$srv_nmfs_trwl_lgth_indicator, rep(0, length.out=extra_years))

    new_data$srv_dom_ll_lgth_indicator <- srv_ll_cal_indic
    new_data$srv_nmfs_trwl_lgth_indicator <- srv_tw_cal_indic

    srv_jpll_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_lgth_indicator, n=extra_years)
    jp_ll_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_lgth_indicator, n=extra_years)

    new_data$srv_jap_ll_lgth_indicator <- srv_jpll_cal_indic
    new_data$srv_jap_fishery_ll_lgth_indicator <- jp_ll_cal_indic

    # Survey new_data
    srv_ll_rpn_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_bio_indicator, n=extra_years)
    srv_ll_rpn_obs <- survey_indices$ll_rpn
    srv_ll_rpn_ses <- model_options$obs_pars$surv_ll$rpn_cv*survey_indices$ll_rpn

    new_data$srv_dom_ll_bio_indicator <- srv_ll_rpn_indic
    new_data$obs_dom_ll_bio <- c(new_data$obs_dom_ll_bio, as.vector(srv_ll_rpn_obs))
    new_data$se_dom_ll_bio <- c(new_data$se_dom_ll_bio, as.vector(srv_ll_rpn_ses))

    srv_jpll_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_bio_indicator, n=extra_years)
    new_data$srv_jap_ll_bio_indicator <- srv_jpll_indic

    # Be careful here. The trawl survey technically only happens ever other
    # year, but its being simulated occurring in every year.
    srv_tw_rpw_indic <- c(new_data$srv_nmfs_trwl_bio_indicator, rep(1, length.out=extra_years))#c(new_data$srv_nmfs_trwl_bio_indicator, rep(c(0, 1), length.out=extra_years))
    srv_tw_rpw_obs <- survey_indices$tw_rpw
    srv_tw_rpw_ses <- model_options$obs_pars$surv_tw$rpw_cv*survey_indices$tw_rpw

    new_data$srv_nmfs_trwl_bio_indicator <- srv_tw_rpw_indic
    new_data$obs_nmfs_trwl_bio <- c(new_data$obs_nmfs_trwl_bio, as.vector(srv_tw_rpw_obs))
    new_data$se_nmfs_trwl_bio <- c(new_data$se_nmfs_trwl_bio, as.vector(srv_tw_rpw_ses))

    ll_cpue_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_cpue_indicator, n=extra_years)
    jp_ll_rpw_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_bio_indicator, n=extra_years)

    new_data$ll_cpue_indicator <- ll_cpue_indic
    new_data$srv_jap_fishery_ll_bio_indicator <- jp_ll_rpw_indic

    # Parameters

    new_parameters <- readRDS(param_file)

    # Update selectivity parameter start values to those used by the 
    # production assessment
    # rows = blocks, col1 = a50, col2 = delta (selex parameters), dim3 = sex (males then females...)
    # LL fish selex logist
    new_parameters$ln_ll_sel_pars <- array(c(1.9380e+000, 1.4919e+000, 1.0631e+000,
                                            -6.9100e-001, -9.7358e-002, -3.7957e-001,
                                            1.42565337612, 1.2222e+000, 6.5829e-001,
                                            -6.9100e-001 , 5.3936e-001, 8.1499e-001), dim = c(3, 2, 2)) 

    # trawl fish selex gamma
    new_parameters$ln_trwl_sel_pars <- array(c(2.1048e+000, 2.3315e+000,
                                                1.7701e+000, 2.3315e+000), dim = c(1,2,2)) 

    # survey domestic selex logist
    new_parameters$ln_srv_dom_ll_sel_pars <- array(c(9.6856e-001, 5.4886e-001,
                                                    8.8394e-001, 8.8394e-001,
                                                    1.0927e+000, 6.3434e-001,
                                                    1.0046e+000, 1.0046e+000), dim = c(2,2,2))

    # japanese ll fish selex logist
    new_parameters$ln_srv_jap_ll_sel_pars <- array(c(1.42565337612, -6.9100e-001,
                                                    1.42565337612, -6.9100e-001), dim = c(1,2,2))

    # nmfs power function trawl survey
    new_parameters$ln_srv_nmfs_trwl_sel_pars <- array(c(-1.31805743083, -0.393941625864), dim = c(1,1,2))
  
    # Update q parameter starting values to those used 
    # by the production assessment
    new_parameters$ln_srv_dom_ll_q <- 1.8582e+000
    new_parameters$ln_srv_jap_ll_q <- 4.6027e+000
    new_parameters$ln_ll_cpue_q <- rep(-6.5299e+000, 3)
    new_parameters$ln_srv_nmfs_trwl_q <- -1.5314e-001
    new_parameters$ln_srv_jap_fishery_ll_q <- 1.5737e+000
    new_parameters$ln_M <- -2.1788e+000

    # Update parameter start values to those used
    # by the production assessment
    new_parameters$ln_mean_rec <- 3.285214
    new_parameters$ln_ll_F_avg <- -3.023593
    new_parameters$ln_trwl_F_avg <- -4.48993

    ln_M_year_devs <- SpatialSablefishAssessment::extend_vec_last_val(new_parameters$ln_M_year_devs, n=extra_years)
    ln_M_age_devs  <- SpatialSablefishAssessment::extend_vec_last_val(new_parameters$ln_M_age_devs, n=extra_years)
    
    ln_rec_dev <- c(new_parameters$ln_rec_dev, rep(0, extra_years))
    ln_ll_F <- c(new_parameters$ln_ll_F_devs, rep(mean(new_parameters$ln_ll_F_devs), extra_years))
    ln_tw_F <- c(new_parameters$ln_trwl_F_devs, rep(mean(new_parameters$ln_trwl_F_devs), extra_years))

    new_parameters$ln_M_year_devs <- ln_M_year_devs
    new_parameters$ln_M_age_devs <- ln_M_age_devs
    new_parameters$ln_rec_dev <- ln_rec_dev
    new_parameters$ln_ll_F_devs <- ln_ll_F
    new_parameters$ln_trwl_F_devs <- ln_tw_F

    # my_model = TMB::MakeADFun(data = new_data,
    #                            parameters = new_parameters,
    #                            DLL = "SpatialSablefishAssessment_TMBExports")

    # mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))
    # mle_report = my_model$report(mle_optim$par)  
    # assessment_ssb <- SpatialSablefishAssessment::get_SSB(mle_report) %>% filter(Year == max(Year)) %>% pull(SSB)

    saveRDS(new_data, paste0("data/sablefish_em_data_curr_",file_suffix,".RDS"))
    saveRDS(new_parameters, paste0("data/sablefish_em_par_curr_",file_suffix,".RDS"))

    capture.output(valid <- SpatialSablefishAssessment::validate_input_data_and_parameters(new_data, new_parameters))

    if(valid){
        return(afscOM::listN(new_data, new_parameters))
    }else{
        print("Something was wrong with the EM data.")
        return(FALSE)
    }
}

#' Simulate OM Observation Data for TMB Assessment Model
#' 
#' Format data and observations from an `afscOM` operating model
#' and update the model data and parameters required by the 
#' `SpatialSablefishAssessment` TMB assessment model. 
#' 
#' This function requires the existence of 'data/sablefish_em_data_2022.RDS',
#' and 'data/sablefish_em_par_2022.RDS' files, which come packaged with the
#' `SablefishMSE` codebase. 
#'
#' @param nyears the number of years of data that will be passed to the model
#' @param dem_params demographic parameter matrices subsetted to 1 year
#' @param land_caa nyears worth of landed catch-at-age data (dimensins [1, nages, nsexes, nregions, nfleets])
#' @param survey_indices nyears worth of survey indices (LL RPN, LL RPW, and TW RPW)
#' @param fxfish_caa_obs nyears worth of catch-at-age observation from the fixed gear fishery
#' @param twfish_caa_obs nyears worth of catch-at-age observation from the trawl gear fishery
#' @param ll_ac_obs nyears worth of age composition observations frmo the longline survey
#' @param tw_ac_obs nyears worth of age composition observations frmo the trawl survey
#' @param model_options list of model options provided to the OM
#' @param added_year the number of new years of data being added (should usually be 1) 
#' @param file_suffix suffix to append to saved outputs
#'
#' @export format_em_data
#'
#' @example
#'
simulate_em_data <- function(
    nyears, 
    dem_params, 
    land_caa, 
    survey_indices, 
    fxfish_caa_obs, 
    twfish_caa_obs, 
    ll_ac_obs, 
    tw_ac_obs, 
    ll_srv_indic,
    tw_srv_indic,
    model_options, 
    added_years=1, 
    file_suffix=""
){

    extra_years <- nyears-63
    year_idxs <- 1:(63+extra_years)

    data_file <- "data/sablefish_em_data_2022.RDS"
    param_file <- "data/sablefish_em_par_2022.RDS"

    new_data <- readRDS(data_file)

    # extract ages
    ages <- as.numeric(colnames(dem_params$waa))
    new_data$ages <- ages

    # extract years
    years <- as.double(1960:(1960+nyears-1))
    new_data$years <- years[year_idxs]

    # extract maturity
    maturity <- dem_params$mat[,,1,1]
    if(is.null(dim(maturity))){
        maturity <- matrix(maturity, ncol=1)
    }else{
        maturity <- t(maturity)
    }
    new_data$maturity <- maturity[,year_idxs]

    # extract mortality
    mortality <- dem_params$mort[,,1,1]
    if(is.null(dim(mortality))){
        mortality <- matrix(mortality, ncol=1)
    }else{
        mortality <- t(mortality)
    }
    new_data$M <- mortality[,year_idxs]

    # extract WAA
    waa <- dem_params$waa[,,,1, drop=FALSE]
    waa_m <- waa[,,2,]
    waa_f <- waa[,,1,]

    if(is.null(dim(waa_m))){
        waa_m <- matrix(waa_m, ncol=1)
        waa_f <- matrix(waa_f, ncol=1)
    }else{
        waa_m <- t(waa_m)
        waa_f <- t(waa_f)
    }

    new_data$male_mean_weight_by_age <- waa_m[,year_idxs]
    new_data$female_mean_weight_by_age <- waa_f[,year_idxs]

    # age-length transition matrix
    male_al_trans_mat <- extend_3darray_last_dim(new_data$male_age_length_transition, n=extra_years)
    female_al_trans_mat <- extend_3darray_last_dim(new_data$female_age_length_transition, n=extra_years)

    new_data$male_age_length_transition <- male_al_trans_mat
    new_data$female_age_length_transition <- female_al_trans_mat

    # proportion male
    prop_male <- rep(0.5, length(year_idxs))
    new_data$proportion_male <- prop_male
    new_data$proportion_male2 <- prop_male

    # spawning time proportion
    spawn_time <- rep(0, length(year_idxs))
    new_data$spawning_time_proportion <- spawn_time

    # recruitment bias ramp
    new_data$end_yr_rec_est_idx <- nyears-1

    # fishery catch vectors
    #landed_caa <- subset_matrix(land_caa[,,,,,drop=FALSE], 1, 6)
    ll_catch <- apply(land_caa[,,,,1, drop=FALSE], 1, sum, na.rm=TRUE)
    tw_catch <- apply(land_caa[,,,,2, drop=FALSE], 1, sum, na.rm=TRUE)

    new_data$ll_fishery_catch <- ll_catch
    new_data$trwl_fishery_catch <- tw_catch

    # fishery selectivity
    ll_sel_yr_indic <- extend_vec_last_val(new_data$ll_sel_by_year_indicator, n=extra_years)
    tw_sel_yr_indic <- extend_vec_last_val(new_data$trwl_sel_by_year_indicator, n=extra_years)
    jp_sel_yr_indic <- extend_vec_last_val(new_data$srv_jap_fishery_ll_sel_by_year_indicator, n=extra_years)

    new_data$ll_sel_by_year_indicator <- ll_sel_yr_indic
    new_data$trwl_sel_by_year_indicator <- tw_sel_yr_indic
    new_data$srv_jap_fishery_ll_sel_by_year_indicator <- jp_sel_yr_indic

    # survey selectivity and catchability
    srv_ll_sel_yr_indic <- extend_vec_last_val(new_data$srv_dom_ll_sel_by_year_indicator, n=extra_years)
    srv_tw_sel_yr_indic <- extend_vec_last_val(new_data$srv_nmfs_trwl_sel_by_year_indicator, n=extra_years)
    srv_jp_sel_yr_indic <- extend_vec_last_val(new_data$srv_jap_ll_sel_by_year_indicator, n=extra_years)

    new_data$srv_dom_ll_sel_by_year_indicator <- srv_ll_sel_yr_indic
    new_data$srv_nmfs_trwl_sel_by_year_indicator <- srv_tw_sel_yr_indic
    new_data$srv_jap_ll_sel_by_year_indicator <- srv_jp_sel_yr_indic

    srv_ll_q_yr_indic <- extend_vec_last_val(new_data$srv_dom_ll_q_by_year_indicator, n=extra_years)
    srv_tw_q_yr_indic <- extend_vec_last_val(new_data$srv_nmfs_trwl_q_by_year_indicator, n=extra_years)
    srv_jp_q_yr_indic <- extend_vec_last_val(new_data$srv_jap_ll_q_by_year_indicator, n=extra_years)
    jp_q_yr_indic <- extend_vec_last_val(new_data$srv_jap_fishery_ll_q_by_year_indicator, n=extra_years)

    new_data$srv_dom_ll_q_by_year_indicator <- srv_ll_q_yr_indic
    new_data$srv_nmfs_trwl_q_by_year_indicator <- srv_tw_q_yr_indic
    new_data$srv_jap_ll_q_by_year_indicator <- srv_jp_q_yr_indic
    new_data$srv_jap_fishery_ll_q_by_year_indicator <- jp_q_yr_indic

    ll_cpue_q_yr_indic <- extend_vec_last_val(new_data$ll_cpue_q_by_year_indicator, n=extra_years)
    new_data$ll_cpue_q_by_year_indicator <- ll_cpue_q_yr_indic

    # catch at age new_data
    ll_caa_indic <- extend_vec_last_val(new_data$ll_catchatage_indicator, n=extra_years)
    ll_caa_indic[length(new_data$ll_catchatage_indicator[year_idxs]):(nyears-1)] <- 1
    new_data$ll_catchatage_indicator <- ll_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_ll_caa <- apply(fxfish_caa_obs[as.logical(new_data$ll_catchatage_indicator),,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_ll_catchatage <- apply(obs_ll_caa, 1, as.double)
    
    # new_data$ll_catchatlgth_indicator
    # new_data$obs_ll_catchatlgth_m

    tw_caa_indic <- extend_vec_last_val(new_data$trwl_catchatlgth_indicator, n=extra_years)
    tw_caa_indic[length(new_data$trwl_catchatlgth_indicator[year_idxs]):(nyears-1)] <- 1
    new_data$trwl_catchatage_indicator <- tw_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_tw_caa <- apply(twfish_caa_obs[as.logical(new_data$trwl_catchatage_indicator),,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_trwl_catchatage <- apply(obs_tw_caa, 1, as.double)

    srv_ll_caa_indic <- new_data$srv_dom_ll_age_indicator
    srv_ll_caa_indic[length(srv_ll_caa_indic)] <- 1
    srv_ll_caa_indic <- extend_vec_last_val(srv_ll_caa_indic, n=extra_years)
    srv_ll_caa_indic[min(length(new_data$srv_dom_ll_age_indicator), length(year_idxs)):(nyears-1)] <- ll_srv_indic
    srv_ll_caa_indic[nyears] <- 0
    new_data$srv_dom_ll_age_indicator <- srv_ll_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_srv_ll_caa <- apply(ll_ac_obs[as.logical(new_data$srv_dom_ll_age_indicator),,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_srv_dom_ll_age <- apply(obs_srv_ll_caa, 1, as.double)

    srv_jpll_caa_indic <- extend_vec_last_val(new_data$srv_jap_ll_age_indicator, n=extra_years)
    new_data$srv_jap_ll_age_indicator <- srv_jpll_caa_indic

    srv_tw_caa_indic <- new_data$srv_nmfs_trwl_age_indicator
    srv_tw_caa_indic[length(srv_tw_caa_indic)] <- 1
    srv_tw_caa_indic <- extend_vec_last_val(srv_tw_caa_indic, n=extra_years)
    srv_tw_caa_indic[min(length(new_data$srv_nmfs_trwl_age_indicator), length(year_idxs)):(nyears-1)] <- tw_srv_indic
    srv_tw_caa_indic[nyears] <- 0
    new_data$srv_nmfs_trwl_age_indicator <- srv_tw_caa_indic

    # srv_tw_caa_indic <- extend_vec_last_val(new_data$srv_nmfs_trwl_age_indicator, n=extra_years)
    # srv_tw_caa_indic[length(new_data$srv_nmfs_trwl_age_indicator[year_idxs]):(nyears-1)] <- 1
    # new_data$srv_nmfs_trwl_age_indicator <- srv_tw_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_srv_tw_caa <- apply(tw_ac_obs[as.logical(new_data$srv_nmfs_trwl_age_indicator),,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_srv_nmfs_trwl_age <- apply(obs_srv_tw_caa, 1, as.double)

    # catch at length
    ll_cal_indic <- extend_vec_last_val(new_data$ll_catchatlgth_indicator, n=extra_years)
    tw_cal_indic <- extend_vec_last_val(new_data$trwl_catchatlgth_indicator, n=extra_years)

    new_data$ll_catchatlgth_indicator <- ll_cal_indic
    new_data$trwl_catchatlgth_indicator <- tw_cal_indic

    if(extra_years <= 0){
        new_data$obs_ll_catchatlgth_m <- new_data$obs_ll_catchatlgth_m[,1:sum(ll_cal_indic)]
        new_data$obs_ll_catchatlgth_f <- new_data$obs_ll_catchatlgth_f[,1:sum(ll_cal_indic)]
        new_data$obs_trwl_catchatlgth_m <- new_data$obs_trwl_catchatlgth_m[,1:sum(tw_cal_indic)]
        new_data$obs_trwl_catchatlgth_f <- new_data$obs_trwl_catchatlgth_m[,1:sum(tw_cal_indic)]
    }

    if(extra_years > 0){
        srv_ll_cal_indic <- c(new_data$srv_dom_ll_lgth_indicator, rep(0, length.out=extra_years))
        srv_tw_cal_indic <- c(new_data$srv_nmfs_trwl_lgth_indicator, rep(0, length.out=extra_years))
    }else{
        srv_ll_cal_indic <- new_data$srv_dom_ll_lgth_indicator[year_idxs]
        srv_tw_cal_indic <- new_data$srv_nmfs_trwl_lgth_indicator[year_idxs]
    }

    if(extra_years <= 0){
        new_data$obs_srv_dom_ll_lgth_m <- new_data$obs_srv_dom_ll_lgth_m[,1:sum(srv_ll_cal_indic)]
        new_data$obs_srv_dom_ll_lgth_f <- new_data$obs_srv_dom_ll_lgth_f[,1:sum(srv_ll_cal_indic)]
        new_data$obs_srv_nmfs_trwl_lgth_m <- new_data$obs_srv_nmfs_trwl_lgth_m[,1:sum(srv_tw_cal_indic)]
        new_data$obs_srv_nmfs_trwl_lgth_f<- new_data$obs_srv_nmfs_trwl_lgth_f[,1:sum(srv_tw_cal_indic)]
    }

    new_data$srv_dom_ll_lgth_indicator <- srv_ll_cal_indic
    new_data$srv_nmfs_trwl_lgth_indicator <- srv_tw_cal_indic

    srv_jpll_cal_indic <- extend_vec_last_val(new_data$srv_jap_ll_lgth_indicator, n=extra_years)
    jp_ll_cal_indic <- extend_vec_last_val(new_data$srv_jap_fishery_ll_lgth_indicator, n=extra_years)

    new_data$srv_jap_ll_lgth_indicator <- srv_jpll_cal_indic
    new_data$srv_jap_fishery_ll_lgth_indicator <- jp_ll_cal_indic

    # Survey new_data

    # srv_ll_caa_indic <- new_data$srv_dom_ll_age_indicator
    # srv_ll_caa_indic[length(srv_ll_caa_indic)] <- 1
    # srv_ll_caa_indic <- extend_vec_last_val(srv_ll_caa_indic, n=extra_years)
    # # srv_ll_caa_indic[length(new_data$srv_dom_ll_age_indicator[year_idxs]):(nyears-1)] <- 1
    # srv_ll_caa_indic[min(length(new_data$srv_dom_ll_age_indicator), length(year_idxs)):(nyears-1)] <- ll_srv_indic
    # srv_ll_caa_indic[nyears] <- 0

    srv_ll_rpn_indic <- new_data$srv_dom_ll_bio_indicator
    srv_ll_rpn_indic <- extend_vec_last_val(srv_ll_rpn_indic, n=extra_years)
    srv_ll_rpn_indic[min(length(new_data$srv_dom_ll_bio_indicator), length(year_idxs)):(nyears)] <- ll_srv_indic

    joint_cv_ll <- model_options$obs_pars$rpn_cv[3]*2.28
    # joint_sd_ll <- sqrt(log(model_options$obs_pars$rpn_cv[3]^2 + 1))*5
    # joint_cv_ll <- sqrt(exp(joint_sd_ll^2) - 1)

    srv_ll_rpn_obs <- survey_indices$rpns[which(srv_ll_rpn_indic == 1),1,1,1,1]
    srv_ll_rpn_ses <- joint_cv_ll*survey_indices$rpns[which(srv_ll_rpn_indic == 1),1,1,1,1]
    srv_ll_rpn_ses[srv_ll_rpn_ses == 0] <- joint_cv_ll 

    new_data$srv_dom_ll_bio_indicator <- srv_ll_rpn_indic
    new_data$obs_dom_ll_bio <- as.vector(srv_ll_rpn_obs)
    new_data$se_dom_ll_bio <- as.vector(srv_ll_rpn_ses)

    srv_jpll_indic <- extend_vec_last_val(new_data$srv_jap_ll_bio_indicator, n=extra_years)
    new_data$srv_jap_ll_bio_indicator <- srv_jpll_indic

    # Be careful here. The trawl survey technically only happens ever other
    # year, but its being simulated occurring in every year.
    # if(extra_years > 0){
    #     srv_tw_rpw_indic <- c(new_data$srv_nmfs_trwl_bio_indicator, rep(1, length.out=extra_years))#c(new_data$srv_nmfs_trwl_bio_indicator, rep(c(0, 1), length.out=extra_years))
    # }else{
    #     srv_tw_rpw_indic <- new_data$srv_nmfs_trwl_bio_indicator[year_idxs]
    # }
    srv_tw_rpw_indic <- new_data$srv_nmfs_trwl_bio_indicator
    srv_tw_rpw_indic <- extend_vec_last_val(srv_tw_rpw_indic, n=extra_years)
    srv_tw_rpw_indic[min(length(new_data$srv_nmfs_trwl_bio_indicator), length(year_idxs)):(nyears)] <- tw_srv_indic

    
    joint_cv_tw <- model_options$obs_pars$rpn_cv[4]*2.28
    # joint_sd_tw <- sqrt(log(model_options$obs_pars$rpn_cv[4]^2 + 1))*5
    # joint_cv_tw <- sqrt(exp(joint_sd_ll^2) - 1)

    srv_tw_rpw_obs <- survey_indices$rpws[which(srv_tw_rpw_indic == 1),1,1,1,2]
    srv_tw_rpw_ses <- joint_cv_tw*survey_indices$rpws[which(srv_tw_rpw_indic == 1),1,1,1,2]
    srv_tw_rpw_ses[srv_tw_rpw_ses == 0] <- joint_cv_tw

    new_data$srv_nmfs_trwl_bio_indicator <- srv_tw_rpw_indic
    new_data$obs_nmfs_trwl_bio <- as.vector(srv_tw_rpw_obs)
    new_data$se_nmfs_trwl_bio <- as.vector(srv_tw_rpw_ses)

    ll_cpue_indic <- extend_vec_last_val(new_data$ll_cpue_indicator, n=extra_years)
    jp_ll_rpw_indic <- extend_vec_last_val(new_data$srv_jap_fishery_ll_bio_indicator, n=extra_years)

    new_data$ll_cpue_indicator <- ll_cpue_indic
    new_data$srv_jap_fishery_ll_bio_indicator <- jp_ll_rpw_indic

    new_data$obs_ll_cpue <- new_data$obs_ll_cpue[1:sum(ll_cpue_indic)]
    new_data$se_ll_cpue <- new_data$se_ll_cpue[1:sum(ll_cpue_indic)]

    # Turn off model weights and enable TMB estimation for all parameters 
    new_data$catch_likelihood <- 0
    new_data$ll_catchatage_comp_likelihood <- 0
    new_data$trwl_catchatage_comp_likelihood <- 0
    new_data$srv_dom_ll_age_comp_likelihood <- 0
    new_data$srv_nmfs_trwl_age_comp_likelihood <- 0
    new_data$srv_dom_ll_bio_likelihood <- 0
    new_data$srv_nmfs_trwl_bio_likelihood <- 0
    new_data$ll_catchatlgth_comp_likelihood <- 1
    new_data$trwl_catchatlgth_comp_likelihood <- 1
    new_data$srv_dom_ll_lgth_comp_likelihood <- 1
    new_data$srv_nmfs_trwl_lgth_comp_likelihood <- 1

    new_data$loglik_wgt_ll_catch <- 10
    new_data$loglik_wgt_trwl_catch <- 10
    new_data$loglik_wgt_ll_catchatage <- 1
    new_data$loglik_wgt_ll_catchatlgth_m <- 1
    new_data$loglik_wgt_ll_catchatlgth_f <- 1
    new_data$loglik_wgt_trwl_catchatlgth_m <- 1
    new_data$loglik_wgt_trwl_catchatlgth_f <- 1
    new_data$loglik_wgt_srv_dom_ll_age <- 1
    new_data$loglik_wgt_srv_dom_ll_lgth_m <- 1
    new_data$loglik_wgt_srv_dom_ll_lgth_f <- 1
    new_data$loglik_wgt_srv_nmfs_trwl_age <- 1
    new_data$loglik_wgt_srv_nmfs_trwl_lgth_m <- 1
    new_data$loglik_wgt_srv_nmfs_trwl_lgth_f <- 1
    new_data$loglik_wgt_trwl_catchatage <- 1
    new_data$loglik_wgt_srv_dom_ll_bio <- 1
    new_data$loglik_wgt_srv_nmfs_trwl_bio <- 1

    # Turn on priors for Q and M
    new_data$mu_srv_dom_ll_q <- 6.412379        # This comes from the 2023 Sablefish assessment
    new_data$sd_srv_dom_ll_q <- 0.01
    new_data$mu_srv_nmfs_trwl_q <- 0.8580096    # This comes from the 2023 Sablefish assessment
    new_data$sd_srv_nmfs_trwl_q <- 0.01
    new_data$loglik_wgt_q_priors <- 1

    new_data$mu_M <- dem_params$mort[1,1,1,1]
    new_data$sd_M <- 8.5481e-002                # This comes from the 2023 Sablefish assessment

    # Parameters

    new_parameters <- readRDS(param_file)

    # Update selectivity parameter start values to those used by the 
    # production assessment
    # rows = blocks, col1 = a50, col2 = delta (selex parameters), dim3 = sex (males then females...)
    # LL fish selex logist
    new_parameters$ln_ll_sel_pars <- array(c(1.9380e+000, 1.4919e+000, 1.0631e+000,
                                            -6.9100e-001, -9.7358e-002, -3.7957e-001,
                                            1.42565337612, 1.2222e+000, 6.5829e-001,
                                            -6.9100e-001 , 5.3936e-001, 8.1499e-001), dim = c(3, 2, 2)) 

    # trawl fish selex gamma
    new_parameters$ln_trwl_sel_pars <- array(c(2.1048e+000, 2.3315e+000,
                                                1.7701e+000, 2.3315e+000), dim = c(1,2,2)) 

    # survey domestic selex logist
    new_parameters$ln_srv_dom_ll_sel_pars <- array(c(9.6856e-001, 5.4886e-001,
                                                    8.8394e-001, 8.8394e-001,
                                                    1.0927e+000, 6.3434e-001,
                                                    1.0046e+000, 1.0046e+000), dim = c(2,2,2))

    # japanese ll fish selex logist
    new_parameters$ln_srv_jap_ll_sel_pars <- array(c(1.42565337612, -6.9100e-001,
                                                    1.42565337612, -6.9100e-001), dim = c(1,2,2))

    # nmfs power function trawl survey
    new_parameters$ln_srv_nmfs_trwl_sel_pars <- array(c(-1.31805743083, -0.393941625864), dim = c(1,1,2))

    # Update q parameter starting values to those used 
    # by the production assessment
    new_parameters$ln_srv_dom_ll_q <- 1.8582e+000
    new_parameters$ln_srv_jap_ll_q <- 4.6027e+000
    new_parameters$ln_ll_cpue_q <- rep(-6.5299e+000, 3)
    new_parameters$ln_srv_nmfs_trwl_q <- -1.5314e-001
    new_parameters$ln_srv_jap_fishery_ll_q <- 1.5737e+000
    new_parameters$ln_M <- -2.1788e+000
    new_parameters$ln_rec_sex_ratio <- log(0.50)

    # Update parameter start values to those used
    # by the production assessment
    new_parameters$ln_mean_rec <- 3.285214
    new_parameters$ln_ll_F_avg <- -3.023593
    new_parameters$ln_trwl_F_avg <- -4.48993

    ln_M_year_devs <- extend_vec_last_val(new_parameters$ln_M_year_devs, n=extra_years)
    ln_M_age_devs  <- extend_vec_last_val(new_parameters$ln_M_age_devs, n=extra_years)

    if(extra_years > 0){
        ln_rec_dev <- c(new_parameters$ln_rec_dev, rep(0, extra_years))
        ln_ll_F <- c(new_parameters$ln_ll_F_devs, rep(mean(new_parameters$ln_ll_F_devs), extra_years))
        ln_tw_F <- c(new_parameters$ln_trwl_F_devs, rep(mean(new_parameters$ln_trwl_F_devs), extra_years))
    }else{
        ln_rec_dev <- new_parameters$ln_rec_dev[year_idxs]
        ln_ll_F <- new_parameters$ln_ll_F_devs[year_idxs]
        ln_tw_F <- new_parameters$ln_trwl_F_devs[year_idxs]
    }
    
    new_parameters$ln_M_year_devs <- ln_M_year_devs
    new_parameters$ln_M_age_devs <- ln_M_age_devs
    new_parameters$ln_rec_dev <- ln_rec_dev
    new_parameters$ln_ll_F_devs <- ln_ll_F
    new_parameters$ln_trwl_F_devs <- ln_tw_F

    capture.output(valid <- SpatialSablefishAssessment::validate_input_data_and_parameters(new_data, new_parameters))

    if(valid){
        return(afscOM::listN(new_data, new_parameters))
    }else{
        print("Something was wrong with the EM data.")
        return(FALSE)
    }

}

#' Simulate OM Observation Data for Sex-Disaggregated TMB Assessment Model
#' 
#' Format data and observations from an `afscOM` operating model
#' and update the model data and parameters required by the 
#' `SpatialSablefishAssessment` TMB assessment model. Age composition data
#' is properly formatted for use in a sex-disaggregated "proportions across"
#' approch as is used in the TMB estimation model.
#' 
#' This function requires the existence of 'data/sablefish_em_data_2022.RDS',
#' and 'data/sablefish_em_par_2022.RDS' files, which come packaged with the
#' `SablefishMSE` codebase. 
#'
#' @param nyears the number of years of data that will be passed to the model
#' @param dem_params demographic parameter matrices subsetted to 1 year
#' @param land_caa nyears worth of landed catch-at-age data (dimensins [1, nages, nsexes, nregions, nfleets])
#' @param survey_indices nyears worth of survey indices (LL RPN, LL RPW, and TW RPW)
#' @param fxfish_caa_obs nyears worth of catch-at-age observation from the fixed gear fishery
#' @param twfish_caa_obs nyears worth of catch-at-age observation from the trawl gear fishery
#' @param ll_ac_obs nyears worth of age composition observations frmo the longline survey
#' @param tw_ac_obs nyears worth of age composition observations frmo the trawl survey
#' @param model_options list of model options provided to the OM
#' @param added_year the number of new years of data being added (should usually be 1) 
#' @param file_suffix suffix to append to saved outputs
#'
#' @export format_em_data
#'
#' @example
#'
simulate_em_data_sex_disaggregate <- function(
    nyears, 
    dem_params, 
    land_caa, 
    survey_indices, 
    fxfish_caa_obs, 
    twfish_caa_obs, 
    ll_ac_obs, 
    tw_ac_obs, 
    ll_srv_indic,
    tw_srv_indic,
    model_options, 
    added_years=1, 
    file_suffix=""
){
    out <- simulate_em_data(nyears, dem_params, land_caa, survey_indices, fxfish_caa_obs, twfish_caa_obs, ll_ac_obs, tw_ac_obs, ll_srv_indic, tw_srv_indic, model_options, added_years, file_suffix)

    out$new_data$obs_ll_catchatage <- t(abind::abind(fxfish_caa_obs[,,2,], fxfish_caa_obs[,,1,], along=2))
    out$new_data$obs_srv_dom_ll_age <- t(abind::abind(ll_ac_obs[,,2,], ll_ac_obs[,,1,], along=2))
    out$new_data$obs_srv_nmfs_trwl_age <- t(abind::abind(tw_ac_obs[,,2,], tw_ac_obs[,,1,], along=2))
    out$new_data$obs_trwl_catchatage <- t(abind::abind(twfish_caa_obs[,,2,], twfish_caa_obs[,,1,], along=2))
    
    # Need to do this for the selectivity estimation to not be weird
    out$new_data$ages <- as.double(1:30)

    # Not using length comps, so simulate catch at age in every year 
    out$new_data$ll_catchatage_indicator <- rep(1, ncol(out$new_data$obs_ll_catchatage))
    out$new_data$obs_ll_catchatage <- t(apply(out$new_data$obs_ll_catchatage, 1, as.double))
    out$new_data$obs_ll_catchatage <- out$new_data$obs_ll_catchatage[,as.logical(out$new_data$ll_catchatage_indicator)]

    out$new_data$obs_srv_dom_ll_age <- t(apply(out$new_data$obs_srv_dom_ll_age, 1, as.double))
    out$new_data$obs_srv_dom_ll_age <- out$new_data$obs_srv_dom_ll_age[,as.logical(out$new_data$srv_dom_ll_age_indicator)]

    out$new_data$obs_srv_nmfs_trwl_age <- t(apply(out$new_data$obs_srv_nmfs_trwl_age, 1, as.double))
    out$new_data$obs_srv_nmfs_trwl_age <- out$new_data$obs_srv_nmfs_trwl_age[,as.logical(out$new_data$srv_nmfs_trwl_age_indicator)]

    # Not using length comps, so simulate catch at age in every year
    out$new_data$trwl_catchatage_indicator <- rep(1, ncol(out$new_data$obs_trwl_catchatage))
    out$new_data$obs_trwl_catchatage <- t(apply(out$new_data$obs_trwl_catchatage, 1, as.double))
    out$new_data$obs_trwl_catchatage <- out$new_data$obs_trwl_catchatage[,as.logical(out$new_data$trwl_catchatage_indicator)]

    out$new_data$trwl_catchatage_covar_structure <- 0

    # Turn off catch at length datasets
    out$new_data$ll_catchatlgth_indicator <- rep(0, length(out$new_data$ll_catchatlgth_indicator))
    out$new_data$trwl_catchatlgth_indicator <- rep(0, length(out$new_data$trwl_catchatlgth_indicator))
    out$new_data$srv_dom_ll_lgth_indicator <- rep(0, length(out$new_data$srv_dom_ll_lgth_indicator))
    out$new_data$srv_nmfs_trwl_lgth_indicator <- rep(0, length(out$new_data$srv_nmfs_trwl_lgth_indicator))


    return(out)
}
