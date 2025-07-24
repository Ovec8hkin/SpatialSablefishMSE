#' Generate Data, Parameter, and Mapping Inputs for SPoRK Estimation Model
#' 
#' @param nyears number of years data observations
#' @param dem_params list of demographic parameters
#' @param agg_land_caa aggregated landed catch at age array [nyears, nages, nsexes, 1, nfleets]
#' @param aggregated_survey_obs list of aggregated survey observations (indices and age compositions)
#' @param model_options OM model options
#' @param r_ISS input sample sizes for age composition data [nyears, nfleets+nsurveys]
#' 
#' @returns list of data, parameters, and parameter mappings for SPoRK TMB model
#' 
#' @export generate_RTMB_inputs
#' 
generate_RTMB_inputs <- function(nyears, dem_params, agg_land_caa, aggregated_survey_obs, model_options, r_ISS){

    # Model Dimensions
    input_list <- SPoRC::Setup_Mod_Dim(
        years=1:nyears,
        ages=1:30,
        lens=seq(41,99,2),
        n_regions=1,
        n_sexes=2,
        n_fish_fleets=2,
        n_srv_fleets = 2
    )

    # Recruitment
    input_list <- SPoRC::Setup_Mod_Rec(
        input_list = input_list,
        do_rec_bias_ramp = 1,
        bias_year = c(length(1960:1979), length(1960:1989), (length(1960:(1960+nyears-1)) - 5), length(1960:(1960+nyears)) - 2) + 1,
        sigmaR_switch = as.integer(length(1960:1975)), 
        dont_est_recdev_last = 0, 
        ln_sigmaR = log(c(0.4, 1.2)),
        rec_model = "mean_rec", # recruitment model
        sigmaR_spec = "fix", 
        sexratio = as.vector(c(0.5, 0.5)), 
        init_age_strc = 1,
        init_F_prop = 0
    )

    # Biologicals
    WAA <- om_to_sporc(dem_params$waa[1:nyears,,,1,drop=FALSE])
    MatAA <- om_to_sporc(dem_params$mat[1:nyears,,,1,drop=FALSE])
    AgeingError <- as.matrix(SPoRC::sgl_rg_sable_data$age_error)
    SizeAgeTrans = SPoRC::sgl_rg_sable_data$SizeAgeTrans
    M <- dem_params$mort[1,1,1,1]
    input_list <- SPoRC::Setup_Mod_Biologicals(
        input_list = input_list,
        WAA = WAA,
        MatAA = MatAA,
        AgeingError = AgeingError,
        SizeAgeTrans = SizeAgeTrans,
        fit_lengths = 0,
        Use_M_prior = 1,
        M_prior = c(0.1, 0.1),
        M_spec = "est_ln_M_only",
        ln_M = log(M),
        M_offset = 0
    )

    # Movement and Tagging
    input_list <- SPoRC::Setup_Mod_Movement(
        input_list = input_list,
        use_fixed_movement = 1,
        Fixed_Movement = NA,
        do_recruits_move = 0
    )

    input_list <-  SPoRC::Setup_Mod_Tagging(
        input_list = input_list,
        UseTagging = 0
    )

    # Catch and Fishing Mortality
    ObsCatch <- array(
        apply(agg_land_caa, c(1, 5), sum), 
        dim=c(input_list$data$n_regions, nyears, input_list$data$n_fish_fleets)
    )
    ObsCatch[,1:3,2] <- NA
    Catch_Type <- array(1, dim=c(length(input_list$data$years), input_list$data$n_fish_fleets))
    UseCatch <- array(1, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
    UseCatch[,1:3,2] <- 0 # Dont fit to 0 catches for trawl fleet in first 3 years

    input_list <- SPoRC::Setup_Mod_Catch_and_F(
        input_list = input_list,
        # Data inputs
        ObsCatch = ObsCatch,
        Catch_Type = Catch_Type,
        UseCatch = UseCatch,
        # Model options
        Use_F_pen = 1,
        sigmaC_spec = 'fix',
        Catch_Constant = c(0, 0)#c(0.01, 0.8)
    )

    ### Fishery Index and Composition Data ---------------------

    # Fishery Index (e.g., none)
    ObsFishIdx <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
    ObsFishIdx_SE <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
    UseFishIdx <- array(0, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))

    # Fishery Age Compositions
    # May need to mess with sex partitioning later
    ObsFishAgeComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    ObsFishAgeComps <- aperm(aggregated_survey_obs$acs[,,,,1:2,drop=FALSE], perm=c(4, 1, 2, 3, 5))
    UseFishAgeComps <- array(1, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
    UseFishAgeComps[1,1:3,2] <- 0
    # Currently assuming equal sample sizes by sex. Fix later.
    ISS_FishAgeComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    ISS_FishAgeComps <- aperm(array(r_ISS[,1:2], dim=c(1, nyears, 2, 2)), c(1, 2, 4, 3)) 
    # aperm(array(model_options$obs_pars$ac_samps[1:2], dim=c(2, nyears, 2, 1)), c(4, 2, 1, 3))
    # aperm(apply(aggregated_survey_obs$acs[,,,,1:2,drop=FALSE], c(1, 3, 4, 5), sum), c(1, )

    ISS_FishAgeComps[,,2,] <- NA

    # Fishery Length Compositions (e.g., none)
    ObsFishLenComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    UseFishLenComps <- array(0, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_fish_fleets))
    ISS_FishLenComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))

    input_list <- SPoRC::Setup_Mod_FishIdx_and_Comps(
        input_list = input_list,
        ObsFishIdx = ObsFishIdx,
        ObsFishIdx_SE = ObsFishIdx_SE,
        UseFishIdx = UseFishIdx,
        ObsFishAgeComps = ObsFishAgeComps,
        UseFishAgeComps = UseFishAgeComps,
        ISS_FishAgeComps = ISS_FishAgeComps,
        ObsFishLenComps = ObsFishLenComps,
        UseFishLenComps = UseFishLenComps,
        ISS_FishLenComps = ISS_FishLenComps,
        fish_idx_type = c("none", "none"),
        FishAgeComps_LikeType = c("Multinomial", "Multinomial"),
        FishLenComps_LikeType = c("none", "none"),
        FishAgeComps_Type = c(paste0("spltRjntS_Year_1-",nyears,"_Fleet_1"), paste0("spltRjntS_Year_1-",nyears,"_Fleet_2")),
        FishLenComps_Type = c(paste0("none_Year_1-",nyears,"_Fleet_1"), paste0("none_Year_1-",nyears,"_Fleet_2")),
    )

    ### Survey Index and Composition Data ---------------------

    # Survey Index (e.g., none)
    ObsSrvIdx <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
    ObsSrvIdx[,,1] <- aggregated_survey_obs$rpns[,1,1,1,1] # Longline survey
    ObsSrvIdx[,,2] <- aggregated_survey_obs$rpws[,1,1,1,2] # Trawl survey

    ObsSrvIdx_SE <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
    #n=10
    # rpn1_idx <- apply(sapply(ObsSrvIdx[,,1], \(x) afscOM::simulate_lognormal_obs(x, model_options$obs_pars$rpn_cv[3], n=n)), 2, sd)#/sqrt(n)
    # rpn2_idx <- apply(sapply(ObsSrvIdx[,,2], \(x) afscOM::simulate_lognormal_obs(x, model_options$obs_pars$rpn_cv[4], n=n)), 2, sd)#/sqrt(n)
    ObsSrvIdx_SE[,,1] <- sqrt(log(model_options$obs_pars$rpn_cv[3]^2+1))
    ObsSrvIdx_SE[,,2] <- sqrt(log(model_options$obs_pars$rpn_cv[4]^2+1))

    UseSrvIdx <- array(1, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))

    # Srvery Age Compositions
    # May need to mess with sex partitioning later
    ObsSrvAgeComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$ages), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    ObsSrvAgeComps <- aperm(aggregated_survey_obs$acs[,,,,3:4,drop=FALSE], perm=c(4, 1, 2, 3, 5))
    UseSrvAgeComps <- array(1, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
    # Currently assuming equal sample sizes by sex. Fix later.
    ISS_SrvAgeComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    ISS_SrvAgeComps <- aperm(array(r_ISS[,3:4], dim=c(1, nyears, 2, 2)), c(1, 2, 4, 3)) 
    ISS_SrvAgeComps[,,2,] <- NA

    # Srvery Length Compositions (e.g., none)
    ObsSrvLenComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), length(input_list$data$lens), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    UseSrvLenComps <- array(0, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_srv_fleets))
    ISS_SrvLenComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))


    input_list <- SPoRC::Setup_Mod_SrvIdx_and_Comps(
        input_list = input_list,
        ObsSrvIdx = ObsSrvIdx,
        ObsSrvIdx_SE = ObsSrvIdx_SE,
        UseSrvIdx = UseSrvIdx,
        ObsSrvAgeComps = ObsSrvAgeComps,
        UseSrvAgeComps = UseSrvAgeComps,
        ISS_SrvAgeComps = ISS_SrvAgeComps,
        ObsSrvLenComps = ObsSrvLenComps,
        UseSrvLenComps = UseSrvLenComps,
        ISS_SrvLenComps = ISS_SrvLenComps,
        srv_idx_type = c("abd", "biom"),
        SrvAgeComps_LikeType = c("Multinomial", "Multinomial"),
        SrvLenComps_LikeType = c("none", "none"),
        SrvAgeComps_Type = c(
            paste0("spltRjntS_Year_1-",nyears,"_Fleet_1"), 
            paste0("spltRjntS_Year_1-",nyears,"_Fleet_2")
        ),
        SrvLenComps_Type = c(
            paste0("none_Year_1-",nyears,"_Fleet_1"), 
            paste0("none_Year_1-",nyears,"_Fleet_2")
        )
    )

    ### Selectivity and Catchability ------------------------

    fsh_sel_blocks_base <- c(
        paste0("Block_1_Year_1-",min(56,nyears),"_Fleet_1"), 
        "",
        "none_Fleet_2"
    )
    if(nyears >= 57){
        fsh_sel_blocks_base[2] <- paste0("Block_2_Year_",min(57,nyears),"-",nyears,"_Fleet_1")
    }else{
        fsh_sel_blocks_base <- fsh_sel_blocks_base[-c(2)]
    }

    input_list <- SPoRC::Setup_Mod_Fishsel_and_Q(
        input_list = input_list,                               
        # Model options
        cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
        # fishery selectivity blocks
        # Fixed gear block 1 1-56, block2 57+; Trawl gear no block
        fish_sel_blocks = fsh_sel_blocks_base, 
        # fishery selectivity form
        # Fixed gear logistic selectivity (k, a50); Trawl gear gamma
        fish_sel_model = c("logist1_Fleet_1", "gamma_Fleet_2"),
        # fishery catchability blocks
        # None as no fishery indices
        fish_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
        # whether to estimate all fixed effects 
        # for fishery selectivity and later modify
        # to fix and share parameters 
        fish_fixed_sel_pars = c("est_all", "est_all"),
        # whether to estimate all fixed effects 
        # for fishery catchability
        fish_q_spec = c("fix", "fix") 
    )

    input_list <- SPoRC::Setup_Mod_Srvsel_and_Q(
        input_list = input_list,
        # Model options
        # survey selectivity, whether continuous time-varying
        cont_tv_srv_sel = c("none_Fleet_1", "none_Fleet_2"),
        # survey selectivity blocks
        srv_sel_blocks = c(
            # paste0("Block_1_Year_1-",min(57,nyears),"_Fleet_1"), 
            # paste0("Block_2_Year_",min(58,nyears),"-",nyears,"_Fleet_1"),
            "none_Fleet_1",
            "none_Fleet_2"
        ),
        # survey selectivity form
        srv_sel_model = c("logist1_Fleet_1","exponential_Fleet_2"),
        # survey catchability blocks
        srv_q_blocks = c("none_Fleet_1", "none_Fleet_2"),
        # whether to estimate all fixed effects 
        # for survey selectivity and later
        # modify to fix/share parameters
        srv_fixed_sel_pars_spec = c("est_all", "est_all"),
        # whether to estimate all 
        # fixed effects for survey catchability
        srv_q_spec = c("est_all", "est_all")
    )
    input_list$par$ln_srv_fixed_sel_pars[,,,2,2] <- log(0.26)

    ### Model Weighting -----------------------------------
    Wt_FishAgeComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    Wt_FishAgeComps[,,1,1] <- 1 #1.018256
    Wt_FishAgeComps[,,1,2] <- 1 #0.6406668

    Wt_FishLenComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_fish_fleets))
    Wt_SrvAgeComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))
    Wt_SrvAgeComps[,,1,1] <- 1 #1.072216
    Wt_SrvAgeComps[,,1,2] <- 1 #1.749187

    Wt_SrvLenComps <- array(NA, dim=c(input_list$data$n_regions, length(input_list$data$years), input_list$data$n_sexes, input_list$data$n_srv_fleets))

    input_list <- SPoRC::Setup_Mod_Weighting(
        input_list = input_list,
        sablefish_ADMB = 0,
        likelihoods = 1,
        Wt_Catch = 1, #50,
        Wt_FishIdx = 0,
        Wt_SrvIdx = 1, #5,
        Wt_Rec = 1, #1.5,
        Wt_F = 1, #0.1,
        Wt_FishAgeComps = Wt_FishAgeComps,
        Wt_FishLenComps = Wt_FishLenComps,
        Wt_SrvAgeComps = Wt_SrvAgeComps,
        Wt_SrvLenComps = Wt_SrvLenComps
    )

    ### Extra Parameter Setting ---------------------------
    input_list$par$ln_sigmaC <- array(log(c(0.02, 0.02)), dim=c(1, 2))

    ### Mapping -------------------------------------------
    # input_list$map$ln_fish_fixed_sel_pars <- factor(c(1:7, 2, rep(c(8,9),2), rep(c(10,9),2)))
    # input_list$map$ln_srv_fixed_sel_pars <-  factor(c(1:3, 2, 4:6, 5,rep(7,4), rep(8, 4)))
    # input_list$map$ln_fish_fixed_sel_pars <- factor(c(1:length(input_list$map$ln_fish_fixed_sel_pars)))
    # input_list$map$ln_srv_fixed_sel_pars <- factor(c(1:length(input_list$map$ln_srv_fixed_sel_pars)))

    return(input_list)

}

om_to_sporc <- function(x){
    return(aperm(x, perm=c(4, 1, 2, 3)))
}
