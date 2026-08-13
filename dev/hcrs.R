tier3 <- function(ref_pts, naa, dem_params, avgrec, cutoff_age=1){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- cutoff_age-1
    ssb <- apply(naa[,a:nages,1,,drop=FALSE]*dem_params$waa[,a:nages,1,1,drop=FALSE]*dem_params$mat[,a:nages,1,1,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$Bref, ref_pts$Fref)
    )
}

chr <- function(ref_pts, naa, dem_params, avgrec){
    return(constant_F(ref_pts$Fref))
}

dynamic_b0 <- function(ref_pts, naa, dem_params, avgrec, cutoff_age=1){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- cutoff_age-1
    ssb <- apply(naa[,a:nages,1,,drop=FALSE]*dem_params$waa[,a:nages,1,1,drop=FALSE]*dem_params$mat[,a:nages,1,1,drop=FALSE], 1, sum)
    
    dep <- ssb/ref_pts$dynb0
    return(
       threshold_f(dep, f_min=0, f_max=ref_pts$Fref, lrp=0.05, urp=0.40)
    )
}


# Going to start an MSE Options list distinct from everything else
mp_base <- setup_mp_options() # get default values
mp_base$management$abc_regflt_apportionment <- abc_regionfleet_allocation
mp_base$management$abc_tac_regflt_reduction <- array(1, dim=c(5, 2), dimnames=dimnames(abc_regionfleet_allocation))
mp_base$management$regflt_tac_utilization <- average_regflt_tac_utilization
mp_base$apportionment$func <- rpw_moving_average
mp_base$apportionment$pars <- list(window_size=5)

#'
#' Alternative Reference Point HCRs
#' 
mp_f40 <- mp_base
mp_f40$name <- "F40"
mp_f40$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_f50 <- mp_base
mp_f50$name <- "F50"
mp_f50$ref_points$spr_target <- c(0.5, 0.5)
mp_f50$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_b30f50 <- mp_base
mp_b30f50$name <- "B30/F50"
mp_b30f50$ref_points$spr_target <- c(0.5, 0.3)
mp_b30f50$hcr <- list(
    func = tier3,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)


#'
#' Stability Constraint HCRs
#'
mp_5perc <- mp_base
mp_5perc$name <- "F40 Global 5% Stability"
mp_5perc$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = array(0.05, dim=c(1, 2)),
        harvest_cap = NA
    ),
    units = "F"
)

mp_10perc_regional <- mp_base
mp_10perc_regional$name <- "F40 Regional 10% Stability"
mp_10perc_regional$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = array(rep(0.10, 5), dim=c(5, 2)), # 10% up and down for each region
        harvest_cap = NA
    ),
    units = "F"
)

mp_10perc_regflt <- mp_base
mp_10perc_regflt$name <- "F40 Fleet 10% Stability"
mp_10perc_regflt$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        # 10% up and down for Fixed Gear in each region, no constraint for Trawl
        max_stability = aperm(array(matrix(c(0.1, 1.0)), dim=c(2, 2, 5)), c(3, 1, 2)), 
        harvest_cap = NA
    ),
    units = "F"
)

mp_10perc_up <- mp_base
mp_10perc_up$name <- "F40 +/- 10% Up"
mp_10perc_up$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = array(c(1.0, 0.10), dim=c(1, 2)), # no constraint down, 10% up
        harvest_cap = NA
    ),
    units = "F"
)


#'
#' Harvest Cap HCRs
#'
mp_20cap <- mp_base
mp_20cap$name <- "20k Harvest Cap"
mp_20cap$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 20
    ),
    units = "F"
)

#'
#' Constant Fishing Mortality Rules
#' 
mp_f40chr <- mp_base
mp_f40chr$name <- "Constant F40"
mp_f40chr$hcr <- list(
    func = chr,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)
mp_f40chr$ref_points$spr_target <- c(0.40, 0.001)

mp_f00chr <- mp_base
mp_f00chr$name <- "No Fishing"
mp_f00chr$hcr <- list(
    func = chr,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "TAC"
)
mp_f00chr$ref_points$spr_target <- c(1, 0.001)

#'
#' Other HCRS
#' 
mp_dynamicB0 <- mp_base
mp_dynamicB0$name <- "Dynamic B0"
mp_dynamicB0$use_dynb0 <- TRUE
mp_dynamicB0$ref_points$spr_target <- c(0.40, 0.40)
mp_dynamicB0$hcr <- list(
    func = dynamic_b0,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_f40hybrid <- mp_base
mp_f40hybrid$name <- "F40 Hybrid"
mp_f40hybrid$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = array(c(1.0, 0.10), dim=c(1, 2)),
        harvest_cap = 20
    ),
    units = "F"
)

generate_full_utilization_hcr <- function(hcr_obj){
    hcr <- hcr_obj
    hcr$name <- paste0(hcr$name, " | Full Utilization")
    hcr$management$regflt_tac_utilization <- full_utilization
    return(hcr)
}

mp_f40_fullutil <- generate_full_utilization_hcr(mp_f40)
mp_f50_fullutil <- generate_full_utilization_hcr(mp_f50)
mp_b30f50_fullutil <- generate_full_utilization_hcr(mp_b30f50)
mp_f40hybrid_fullutil <- generate_full_utilization_hcr(mp_f40hybrid)
mp_db0_fullutil <- generate_full_utilization_hcr(mp_dynamicB0)
mp_20cap_fullutil <- generate_full_utilization_hcr(mp_20cap)
mp_f00_fullutil <- generate_full_utilization_hcr(mp_f00chr)
