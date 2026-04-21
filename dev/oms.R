# Define recruitment to occur via historical resampling
sable_om <- readRDS("data/spatial_sablefish_om.RDS")
sable_om$model_options$obs_pars$catch_cv <- c(0.02, 0.02, 0, 0)

assessment <- dget(file.path(here::here(), "data/spatial_sablefsh_inputs.rdat"))
hist_recruits <- t(assessment$recruitment)

dp_y <- afscOM::subset_dem_params(sable_om$dem_params, 64, d=1, drop=FALSE)
joint_selret <- calculate_joint_selret(
    sel = dp_y$sel,
    ret = dp_y$ret,
    prop_fs = c(0.80, 0.20)
)

sbpr <- compute_sbpr(
    nages=30,
    mort = dp_y$mort[,,1,1],
    mat = dp_y$mat[,,1,1],
    waa = dp_y$waa[,,1,1],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    F = 0
)

om_base <- sable_om
om_base$dem_params$sel[1:56,,,,1] <- om_base$dem_params$sel[57:112,,,,1] # Remove selectivity block


### Alternative Movement Scenarios

om_climate_movement <- om_base
om_climate_movement$dem_params$movement[,,55:200,,] <- create_tv_movement(time_trend=30, nyears=length(55:200))


# Beverton-Holt recruitment
om_bh_recruit <- om_base
om_bh_recruit$name <- "Beverton-Holt Recruitment"
om_bh_recruit$recruitment$func <- beverton_holt
om_bh_recruit$recruitment$pars <- list(
    h = 0.85, # could do 0.80
    R0 = 18,
    S0 = sbpr*18,
    sigR = 0.70
)
om_bh_recruit$recruitment$apportionment$func <- proximity_resample_recruit_apportionment
om_bh_recruit$recruitment$apportionment$pars <- list(
    hist_recruits = hist_recruits
)

# Beverton-Holt Cyclic Recruitment
om_bhcyclic_recruit <- om_base
om_bhcyclic_recruit$name <- "Beverton-Holt Cyclic Recruitment"
om_bhcyclic_recruit$recruitment$func <- bevholt_regimes
om_bhcyclic_recruit$recruitment$pars <- list(
    h = 0.85,
    sbpr = sbpr,
    R0 = c(70, 10),
    sigR = c(0.70, 0.70),
    nyears = 500,
    regime_length = c(5, 20),
    starting_regime = 0
)
om_bhcyclic_recruit$recruitment$apportionment$func <- proximity_resample_recruit_apportionment
om_bhcyclic_recruit$recruitment$apportionment$pars <- list(
    hist_recruits = hist_recruits
)

om_crash <- om_base
om_crash$name <- "Crash Recruitment"
om_crash$recruitment$func <- bevholt_regimes
om_crash$recruitment$pars <- list(
    h = 0.85,
    sbpr = sbpr,
    R0 = c(18, 3.5),
    sigR = c(0.70, 0.70),
    nyears = 500,
    regime_length = c(12, 25),
    starting_regime = 0
)
om_crash$recruitment$apportionment$func <- proximity_resample_recruit_apportionment
om_crash$recruitment$apportionment$pars <- list(
    hist_recruits = hist_recruits
)


om_agemove_bhrec <- om_bh_recruit
om_agemove_bhrec$name <- "BH Recruit | AB Move"

om_agemove_regimerec <- om_bhcyclic_recruit
om_agemove_regimerec$name <- "Regime Recruit | AB Move"

om_agemove_crashrec <- om_crash
om_agemove_crashrec$name <- "Crash Recruit | AB Move"

om_climatemove_bhrec <- om_bh_recruit
om_climatemove_bhrec$dem_params$movement <- om_climate_movement$dem_params$movement
om_climatemove_bhrec$name <- "BH Recruit | Climate Move"

om_climatemove_regimerec <- om_bhcyclic_recruit
om_climatemove_regimerec$dem_params$movement <- om_climate_movement$dem_params$movement
om_climatemove_regimerec$name <- "Regime Recruit | Climate Move"

om_climatemove_crashrec <- om_crash
om_climatemove_crashrec$dem_params$movement <- om_climate_movement$dem_params$movement
om_climatemove_crashrec$name <- "Crash Recruit | Climate Move"