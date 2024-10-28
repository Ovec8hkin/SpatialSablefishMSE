
# Define recruitment to occur via historical resampling
assessment <- dget("data/sablefish_assessment_2023.rdat")
hist_recruits <- assessment$natage.female[,1]*2

dp_y <- afscOM::subset_dem_params(sable_om$dem_params, 64, d=1, drop=FALSE)
joint_selret <- calculate_joint_selret(
    sel = dp_y$sel,
    ret = dp_y$ret,
    prop_fs = c(0.80, 0.20)
)
rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = 0.40
)
ref_naa <- compute_naapr(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    F = rp$Fref
)
sbpr <- compute_sbpr(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    F = 0
)



# OM1: Normal recruitment
om_rand_recruit <- sable_om
om_rand_recruit$name <- "Random Recruitment"
om_rand_recruit$recruitment$func <- resample_recruits
om_rand_recruit$recruitment$pars <- list(
    hist_recruits = hist_recruits,
    nyears = 10*nyears
)

# OM2: Cyclic recruiment
om_cylic_recruit <- sable_om
om_cylic_recruit$name <- "Cyclic Recruitment"
om_cylic_recruit$recruitment$func <- resample_regime_recruits
om_cylic_recruit$recruitment$pars <- list(
    regime1_recruits = hist_recruits[hist_recruits <= 35],
    regime2_recruits = hist_recruits[hist_recruits > 35],
    nyears = 10*nyears,
    regime_length = c(20, 5),
    starting_regime = 0
)

# OM3: B-H recruitment
om_bh_recruit <- sable_om
om_bh_recruit$name <- "Beverton-Holt Recruitment"
om_bh_recruit$recruitment$func <- beverton_holt
om_bh_recruit$recruitment$pars <- list(
    h = 0.85,
    R0 = 25,
    S0 = sbpr*25,
    sigR = 1.20
)

# OM4: Low recruitment
om_low_recruit <- sable_om
om_low_recruit$name <- "Low Recruitment"
om_low_recruit$recruitment$func <- regime_recruits
om_low_recruit$recruitment$pars <- list(
    mus = c(mean(hist_recruits[hist_recruits <= 35])),
    cvs = c(sd(hist_recruits[hist_recruits <= 35])/mean(hist_recruits[hist_recruits <= 35])),
    nyears = 10*nyears,
    regime_length = c(nyears*10),
    starting_regime = 0
)

om_high_recruit <- sable_om
om_high_recruit$name <- "High Recruitment"
om_high_recruit$recruitment$func <- regime_recruits
om_high_recruit$recruitment$pars <- list(
    mus = c(mean(hist_recruits[hist_recruits > 35])),
    cvs = c(sd(hist_recruits[hist_recruits > 35])/mean(hist_recruits[hist_recruits > 35])),
    nyears = 10*nyears,
    regime_length = c(nyears*10),
    starting_regime = 0
)

#OM6: Crash recruitment
om_crash_recruit <- sable_om
om_crash_recruit$name <- "Crash Recruitment"
om_crash_recruit$recruitment$func <- recruits_crash
om_crash_recruit$recruitment$pars <- list(
    crash_start_year = (nyears-64)/2,
    crash_length = 20,
    crash_value = min(hist_recruits),
    hist_recruits = hist_recruits,
    nyears = 10*nyears
)

# OM3: B-H recruitment
om_bhcyclic_recruit <- sable_om
om_bhcyclic_recruit$name <- "Beverton-Holt Cyclic Recruitment"
om_bhcyclic_recruit$recruitment$func <- bevholt_regimes
om_bhcyclic_recruit$recruitment$pars <- list(
    h = 0.85,
    sbpr = sbpr,
    R0 = c(12.5, 90),
    sigR = c(1.20, 1.20),
    nyears = 10*nyears,
    regime_length = c(20, 5),
    starting_regime = 0
)