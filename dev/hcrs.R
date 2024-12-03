tier3 <- function(ref_pts, naa, dem_params, avgrec, cutoff_age=1){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- cutoff_age-1
    ssb <- apply(naa[,a:nages,1,]*dem_params$waa[,a:nages,1,,drop=FALSE]*dem_params$mat[,a:nages,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$Bref, ref_pts$Fref)
    )
}

new_hcr <- function(ref_pts, naa, dem_params, avgrec, cutoff_age=1, alpha=0.05, beta=0.34){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- cutoff_age-1
    ssb <- apply(naa[,a:nages,1,]*dem_params$waa[,a:nages,1,,drop=FALSE]*dem_params$mat[,a:nages,1,,drop=FALSE], 1, sum)
    return(
        threshold_cap(ssb, B_ref=ref_pts$Bref, F_ref=ref_pts$Fref, alpha=alpha, beta=beta)
    )
}

pfmc4010 <- function(ref_pts, naa, dem_params, avgrec, pstar=0.45, OFLsigma=0.32){

    joint_selret <- calculate_joint_selret(dem_params$sel, dem_params$ret, c(0.80, 0.20))
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    exploit_bio <- naa*dem_params$waa*joint_selret$sel
    dep <- ssb/ref_pts$B0

    # From Maia
    ABC <- afscOM::F_to_mu(ref_pts$Fref)*sum(exploit_bio)
        # ABC <- ref_pts$Fref
    # ABC <- OFL*exp(qnorm(pstar, 0, OFLsigma))
    if(dep < 0.1){TAC <- 0}
    if(dep > 0.40) {TAC <- ABC}
    if(dep <= 0.40 & dep >= 0.1){
    #   TAC <- ref_pts$Fref*ssb*((ssb-0.1*ref_pts$B0)/(0.4*ref_pts$B0 -0.1*ref_pts$B0))
        TAC <- ABC*(dep-0.1)/(0.4-0.1)
    }

    return(TAC)
}

bc_sable <- function(ref_pts, naa, dem_params, avgrec){
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    dep <- ssb/ref_pts$Bref

    # From Maia
    ABC <- 0.055
    # ABC <- OFL*exp(qnorm(pstar, 0, OFLsigma))
    if(dep < 0.4){TAC <- 0}
    if(dep > 0.60) {TAC <- ABC}
    if(dep <= 0.60 & dep >= 0.40){
    #   TAC <- ref_pts$Fref*ssb*((ssb-0.1*ref_pts$B0)/(0.4*ref_pts$B0 -0.1*ref_pts$B0))
        TAC <- ABC*(dep-0.4)/(0.6-0.4)
    }

    return(TAC)
}

age_structure_percentage <- function(ref_pts, naa, dem_params, avgrec, desired_abi, ref_naa){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- 2-1
    ssb <- apply(naa[,a:nages,1,]*dem_params$waa[,a:nages,1,,drop=FALSE]*dem_params$mat[,a:nages,1,,drop=FALSE], 1, sum)
    f_default <- npfmc_tier3_F(ssb, ref_pts$Bref, ref_pts$Fref)

    abi_y <- abi(naa[,,1,], ref_naa, threshold = 0.90)
    abi_perc <- abi_y/desired_abi

    new_f <- f_default * abi_perc
    new_f <- ifelse(new_f > 0.95*ref_pts$Fmax & new_f > ref_pts$Fref, ref_pts$Fmax*0.95, new_f)
    return(new_f)
}

age_structure_rt_reduction <- function(ref_pts, naa, dem_params, ref_naa, avgrec, breakpoints=c(0.4, 0.8, 1.5), levels=c(0.75, 0.85, 1.0)){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- 2-1
    ssb <- apply(naa[,a:nages,1,]*dem_params$waa[,a:nages,1,,drop=FALSE]*dem_params$mat[,a:nages,1,,drop=FALSE], 1, sum)

    abi_y <- abi(naa[,,1,], ref_naa, threshold = 0.90)

    reduction_factor <- 1
    if(abi_y < breakpoints[1]){
        reduction_factor <- levels[1]
    }else if((abi_y > breakpoints[1] && abi_y < breakpoints[2]) | abi_y > breakpoints[3]){
        reduction_factor <- levels[2]
    }else{
        reduction_factor <- levels[3]
    }


    new_fmax <- ref_pts$Fref * reduction_factor
    new_fmax <- ifelse(new_fmax > ref_pts$Fmax*0.99, ref_pts$Fmax*0.99, new_fmax)

    f <- npfmc_tier3_F(ssb, ref_pts$Bref, new_fmax)

    return(f)

}

ref_point_switch <- function(ref_pts, naa, dem_params, avgrec, abi_switch, ref_naa, new_refpts){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- 2-1
    ssb <- apply(naa[,a:nages,1,]*dem_params$waa[,a:nages,1,,drop=FALSE]*dem_params$mat[,a:nages,1,,drop=FALSE], 1, sum)

    abi_y <- abi(naa[,,1,], ref_naa, threshold = 0.90)
    if(abi_y < abi_switch){
        joint_selret <- calculate_joint_selret(dem_params$sel, dem_params$ret, c(0.80, 0.20))
        ref_pts <- calculate_ref_points(
            nages = nages,
            mort = dem_params$mort[,,1,],
            mat = dem_params$mat[,,1,],
            waa = dem_params$waa[,,1,],
            sel =  joint_selret$sel[,,1,,drop=FALSE],
            ret = joint_selret$ret[,,1,,drop=FALSE],
            avg_rec = avgrec/2,
            spr_target = new_refpts
        )
    }

    f_default <- npfmc_tier3_F(ssb, ref_pts$Bref, ref_pts$Fref)
    return(f_default)
}

chr <- function(ref_pts, naa, dem_params, avgrec){
    return(constant_F(ref_pts$Fref))
}

# Going to start an MSE Options list distinct from everything else
mp_base <- setup_mp_options() # get default values
mp_base$management$tac_land_reduction <- 1

# mp_base$management$tac_land_reduction <- list(
#     func = stairstep_attainment,
#     pars = list(
#         breakpoints = c(20, 30),
#         levels = c(0.874, 0.786, 0.647),
#         phase_ins = 2
#     )
# )

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

mp_b40f55 <- mp_base
mp_b40f55$name <- "B40/F55"
mp_b40f55$ref_points$spr_target <- c(0.55, 0.4)
mp_b40f55$hcr <- list(
    func = tier3,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_b30f40 <- mp_base
mp_b30f40$name <- "B30/F40"
mp_b30f40$ref_points$spr_target <- c(0.4, 0.3)
mp_b30f40$hcr <- list(
    func = tier3,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_b50f40 <- mp_base
mp_b50f40$name <- "B50/F40"
mp_b50f40$ref_points$spr_target <- c(0.4, 0.5)
mp_b50f40$hcr <- list(
    func = tier3,
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

mp_b40f50 <- mp_base
mp_b40f50$name <- "B40/F50"
mp_b40f50$ref_points$spr_target <- c(0.5, 0.4)
mp_b40f50$hcr <- list(
    func = tier3,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_f40cons <- mp_base
mp_f40cons$name <- "F40 Conservative"
mp_f40cons$hcr <- list(
    func = new_hcr,
    extra_pars = NA,
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
mp_5perc$name <- "F40 +/- 5%"
mp_5perc$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = 0.05,
        harvest_cap = NA
    ),
    units = "F"
)

mp_10perc <- mp_base
mp_10perc$name <- "F40 +/- 10%"
mp_10perc$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = 0.10,
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
        max_stability = c(1.0, 0.10), # no constraint down, 10% up
        harvest_cap = NA
    ),
    units = "F"
)

mp_15perc <- mp_base
mp_15perc$name <- "F40 +/- 15%"
mp_15perc$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = 0.15,
        harvest_cap = NA
    ),
    units = "F"
)

mp_25perc <- mp_base
mp_25perc$name <- "F40 +/- 25%"
mp_25perc$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = 0.25,
        harvest_cap = NA
    ),
    units = "F"
)


#'
#' Harvest Cap HCRs
#'
mp_00cap <- mp_base
mp_00cap$name <- "0k Harvest Cap"
mp_00cap$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 0
    ),
    units = "F"
)
mp_00cap$management$tac_land_reduction <- 1 

mp_10cap <- mp_base
mp_10cap$name <- "10k Harvest Cap"
mp_10cap$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 10
    ),
    units = "F"
)
mp_10cap$management$tac_land_reduction <- 1

mp_15cap <- mp_base
mp_15cap$name <- "15k Harvest Cap"
mp_15cap$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 15
    ),
    units = "F"
)
mp_15cap$management$tac_land_reduction <- 1

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
mp_20cap$management$tac_land_reduction <- 1

mp_25cap <- mp_base
mp_25cap$name <- "25k Harvest Cap"
mp_25cap$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 25
    ),
    units = "F"
)
mp_25cap$management$tac_land_reduction <- 1
# mp_25cap$management$tac_land_reduction <- list(
#     func = stairstep_attainment,
#     pars = list(
#         breakpoints = c(20, 30),
#         levels = c(1, 0.874, 0.786),
#         phase_ins = 0
#     )
# )



mp_20cap_cons <- mp_base
mp_20cap_cons$name <- "B30 20k Harvest Cap"
mp_20cap_cons$ref_points$spr_target <- c(0.4, 0.3)
mp_20cap_cons$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 20
    ),
    units = "F"
)
mp_20cap_cons$management$tac_land_reduction <- 1

mp_25cap_cons <- mp_base
mp_25cap_cons$name <- "B30 25k Harvest Cap"
mp_25cap_cons$ref_points$spr_target <- c(0.40, 0.30)
mp_25cap_cons$hcr <- list(
    func = tier3,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = 25
    ),
    units = "F"
)
mp_25cap_cons$management$tac_land_reduction <- 1
# mp_25cap_cons$management$tac_land_reduction <- list(
#     func = stairstep_attainment,
#     pars = list(
#         breakpoints = c(20, 30),
#         levels = c(1, 0.874, 0.786),
#         phase_ins = 0
#     )
# )

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

mp_f50chr <- mp_base
mp_f50chr$name <- "Constant F50"
mp_f50chr$hcr <- list(
    func = chr,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)
mp_f50chr$ref_points$spr_target <- c(0.50, 0.001)


mp_f55chr <- mp_base
mp_f55chr$name <- "Constant F55"
mp_f55chr$hcr <- list(
    func = chr,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)
mp_f55chr$ref_points$spr_target <- c(0.55, 0.001)

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
#' Other HCRs
#' 
mp_agedevs <- mp_base
mp_agedevs$name <- "Age Structure Deviations"
mp_agedevs$hcr <- list(
    func = age_structure_percentage,
    extra_pars = list(
        desired_abi = 1.5,
        ref_naa = ref_naa
    ),
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_ageblocks <- mp_base
mp_ageblocks$name <- "Age Structure Blocks"
mp_ageblocks$hcr <- list(
    func = age_structure_rt_reduction,
    extra_pars = list(
        ref_naa = ref_naa,
        breakpoints=c(0.4, 0.8, 10), 
        levels=c(0.50, 0.75, 1.0)
    ),
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_ageswitch <- mp_base
mp_ageswitch$name <- "Age Structure RP Switch"
mp_ageswitch$hcr <- list(
    func = ref_point_switch,
    extra_pars = list(
        abi_switch = 0.8,
        ref_naa = ref_naa,
        new_refpts = c(0.50, 0.40)
    ),
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)

mp_pfmc4010 <- mp_base
mp_pfmc4010$name <- "PFMC 40-10"
mp_pfmc4010$ref_points$spr_target <- c(0.45, 0.40)
mp_pfmc4010$hcr <- list(
    func = pfmc4010,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "TAC"
)

mp_bcsable <- mp_base
mp_bcsable$name <- "British Columbia"
mp_bcsable$ref_points$spr_target <- c(0.45, 0.45)
mp_bcsable$hcr <- list(
    func = bc_sable,
    extra_pars = NA,
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    units = "F"
)
