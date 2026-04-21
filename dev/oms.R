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

# No movement
om_nomove <- om_base
om_nomove$dem_params$movement[,,55:200,,] <- array(diag(1, nrow=5, ncol=5), dim=dim(om_base$dem_params$movement), dimnames=dimnames(om_base$dem_params$movement))[,,55:200,,]

# Fully mixed
om_panmictic <- om_base
om_panmictic$dem_params$movement[,,55:200,,] <- array(0.20, dim=dim(om_base$dem_params$movement), dimnames=dimnames(om_base$dem_params$movement))[,,55:200,,]

om_climate_movement <- om_base
om_climate_movement$dem_params$movement[,,55:200,,] <- create_tv_movement(time_trend=30, nyears=length(55:200))

# Normal recruitment
om_rand_recruit <- om_base
om_rand_recruit$name <- "Random Recruitment"
om_rand_recruit$recruitment$func <- resample_recruits_spatial
om_rand_recruit$recruitment$pars <- list(
    hist_recruits = hist_recruits,
    nyears = 500
)

# Beverton-Holt recruitment
om_bh_recruit <- om_base
om_bh_recruit$name <- "Beverton-Holt Recruitment"
om_bh_recruit$recruitment$func <- beverton_holt
om_bh_recruit$recruitment$pars <- list(
    h = 0.85, # could do 0.80
    R0 = 15,
    S0 = sbpr*15,
    sigR = 1.20
)
om_bh_recruit$recruitment$apportionment$func <- resample_recruit_apportionment
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
    R0 = c(12.5, 50),
    sigR = c(1.20, 1.20),
    nyears = 500,
    regime_length = c(20, 5),
    starting_regime = 0
)
om_bhcyclic_recruit$recruitment$apportionment$func <- resample_recruit_apportionment
om_bhcyclic_recruit$recruitment$apportionment$pars <- list(
    hist_recruits = hist_recruits
)

# Immediate Recruitment Crash
om_immcrash_recruit <- om_base
om_immcrash_recruit$name <- "Immediate Crash Recruitment"
om_immcrash_recruit$recruitment$func <- recruits_crash
om_immcrash_recruit$recruitment$pars <- list(
    crash_start_year = 1,
    crash_length = 20,
    crash_value = min(hist_recruits),
    hist_recruits = hist_recruits,
    nyears = 500
)
om_immcrash_recruit$recruitment$apportionment$func <- resample_recruit_apportionment
om_immcrash_recruit$recruitment$apportionment$pars <- list(
    hist_recruits = hist_recruits
)


om_crash2 <- om_bhcyclic_recruit
om_crash2$recruitment$pars$h <- c(0.85, 0.85)
om_crash2$recruitment$pars$R0 <- c(15, 3.5)
om_crash2$recruitment$pars$regime_length <- c(25, 30)
om_crash2$name <- "Crash Recruitment"

om_cycle_low <- om_bhcyclic_recruit
om_cycle_low$recruitment$pars$h <- c(0.85,0.85)
om_cycle_low$recruitment$pars$R0 <- c(50, 5.5)
om_cycle_low$recruitment$pars$regime_length <- c(5, 20)
om_cycle_low$name <- "Low Regime Recruitment"
