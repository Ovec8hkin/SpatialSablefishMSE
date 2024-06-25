#' Run Management Strategy Evaluation
#' 
#' Run an MSE simulation loop.
#'
#' @param om path to a saved OM .RDS file (should contain
#' demagraphic parameters list and model options list)
#' @param hcr a function to compute the allowable F in the next year
#' @param ... parameters to pass to the `hcr` function
#' @param nyears_input number of years to simulate forward (will override the dimensions
#' of the demagraphic parameters matrices defined in the OM)
#' @param spinup_years number of years before estimation process should begin
#' @param seed random seed
#'
#' @export run_mse
#'
#' @example
#'
run_mse <- function(om, mp, ..., run_estimation=TRUE, nyears_input=NA, spinup_years=64, seed=1120, file_suffix=""){
   
    assessment <- dget("data/sablefish_assessment_2023.rdat")
   
    # Load OM parameters into global environment
    list2env(om, env=environment())

    # Load OM dimensions into global environment
    list2env(afscOM::get_model_dimensions(dem_params$sel), env=environment())
    nsurveys <- ifelse(model_options$simulate_observations, get_model_dimensions(dem_params$surv_sel)$nfleets, 0)

    if(!is.na(nyears_input)){
        nyears <- nyears_input
    }
    
    dimension_names <- list(
        "time" = 1:nyears,
        "age"  = 2:31,
        "sex"  = c("F", "M"),
        "region" = "alaska",
        "fleet" = c("Fixed", "Trawl")
    )

    landings <- rep(0, nyears)
    hcr_F <- rep(0, nyears)
    out_f <- rep(0, nyears) # vector to store outputted F
    landings[1:64] <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])

    #' 6. Setup empty array to collect derived quantities from the OM
    #' It is left to the user to decide what information to store, and
    #' in what format they would like to store it.
    #'
    #' Here, we are storing all OM outputs as they are returned from the
    #' OM.
    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
    faa_est     = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets), dimnames=dimension_names)
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska"))
    naa_est     = array(NA, dim=c(nyears,   nages, nsexes, nregions), dimnames=list("time"=1:(nyears),   "age"=2:31, "sex"=c("F", "M"), "region"="Alaska"))

    abc         = array(NA, dim=c(nyears+1, 1, 1, 1), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska"))
    tac         = array(NA, dim=c(nyears+1, 1, 1, 1), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska"))
    exp_land    = array(NA, dim=c(nyears+1, 1, 1, 1), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska"))
    hcr_f       = array(NA, dim=c(nyears,   1, 1, 1), dimnames=list("time"=1:nyears,     1, 1, "region"="Alaska"))
    out_f       = array(NA, dim=c(nyears,   1, 1, 1), dimnames=list("time"=1:nyears,     1, 1, "region"="Alaska"))

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets))

    survey_preds <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    survey_obs <- list(
        rpns = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        rpws = array(NA, dim=c(nyears, 1, 1, nregions, nsurveys)),
        acs  = array(NA, dim=c(nyears, nages, nsexes, nregions, nsurveys+nfleets))
    )

    model_outs = list(
        mods = vector("list", length=c(nyears-spinup_years)),
        fits = vector("list", length=c(nyears-spinup_years)),
        reps = vector("list", length=c(nyears-spinup_years))
    )

    naa[1,,,] = init_naa

    set.seed(seed)
    hist_recruitment <- assessment$natage.female[,1]*2
    projected_recruitment <- do.call(recruitment$func, c(recruitment$pars, list(seed=seed)))
    if(!is.function(projected_recruitment)){
        full_recruitment <- c(hist_recruitment, projected_recruitment)
    }else{
        full_recruitment <- rep(NA, nyears)
        full_recruitment[1:length(hist_recruitment)] <- hist_recruitment
    }
    
    for(y in 1:nyears_input){
        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp_y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        removals_input <- landings[y]
        fleet.props <- subset_matrix(model_options$fleet_apportionment, y, d=1, drop=FALSE)
        region_props <- as.matrix(1)
        rec_props <- as.matrix(1)

        # will work for any recruitment function that only requires
        # ssb as a yearly input (beverton holt and ricker should work fine)
        if(is.na(full_recruitment[y+1])){
            ssb <- sum(naa[y,,1,,drop=FALSE]*dp_y$waa[,,1,]*dp_y$mat[,,1,])
            full_recruitment[y+1] <- projected_recruitment(ssb) 
        }

        prev_naa <- naa[y,,,, drop = FALSE]
        out_vars <- project(
            removals = removals_input,
            dem_params=dp_y,
            prev_naa=prev_naa,
            recruitment=full_recruitment[y+1],
            fleet_props = fleet.props,
            region_props = region_props,
            rec_props = rec_props,
            options=model_options
        )

        # update state
        land_caa[y,,,,] <- out_vars$land_caa_tmp
        disc_caa[y,,,,] <- out_vars$disc_caa_tmp
        caa[y,,,,] <- out_vars$caa_tmp
        faa[y,,,,] <- out_vars$faa_tmp
        naa[y+1,,,] <- out_vars$naa_tmp

        f[y,,,,] <- out_vars$F_f_tmp

        survey_preds$rpns[y,,,,] <- out_vars$survey_preds$rpns
        survey_preds$rpws[y,,,,] <- out_vars$survey_preds$rpws
        survey_preds$acs[y,,,,]  <- out_vars$survey_preds$acs

        survey_obs$rpns[y,,,,] <- out_vars$survey_obs$rpns
        survey_obs$rpws[y,,,,] <- out_vars$survey_obs$rpws
        survey_obs$acs[y,,,,]  <- out_vars$survey_obs$acs

        
        if((y+1) > spinup_years){

            naa_proj <- out_vars$naa_tmp
            rec <- full_recruitment[1:y]
            sel <- dp_y$sel
            prop_fs <- apply(out_vars$faa_tmp[,,,1,, drop=FALSE], 5, max)/sum(apply(out_vars$faa_tmp[,,,1,, drop=FALSE], 5, max))

            if(run_estimation){
                # Do all of the data formatting and running
                # of the TMB Sablefish model
                assess_inputs <- simulate_em_data_sex_disaggregate(
                    nyears = y,
                    dem_params = afscOM::subset_dem_params(om$dem_params, 1:y, d=1, drop=FALSE),
                    land_caa = land_caa[1:y,,,,,drop=FALSE],
                    survey_indices = afscOM::subset_dem_params(survey_obs, 1:y, d=1, drop=FALSE),
                    fxfish_caa_obs = afscOM::subset_matrix(survey_obs$acs[1:y,,,,1,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    twfish_caa_obs = afscOM::subset_matrix(survey_obs$acs[1:y,,,,2,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    ll_ac_obs = afscOM::subset_matrix(survey_obs$acs[1:y,,,,3,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    tw_ac_obs = afscOM::subset_matrix(survey_obs$acs[1:y,,,,4,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    model_options = om$model_options,
                    added_years = y-spinup_years,
                    file_suffix = y
                )


                mod_out <- fit_TMB_model(
                    assess_inputs$new_data, 
                    assess_inputs$new_parameters, 
                    model_name = "CurrentAssessmentDisaggregated"
                )  
                mod_report <- mod_out$report
                ssb <- SpatialSablefishAssessment::get_SSB(mod_report) %>% filter(Year == max(Year)) %>% pull(SSB)
                rec <- SpatialSablefishAssessment::get_recruitment(mod_report) %>% filter(Year != max(Year)) %>% pull(Recruitment)
                selex <- SpatialSablefishAssessment::get_selectivities(mod_report)

                # Store assessment estimates of age composition
                # for comparing EM and OM
                if(y == spinup_years){
                    naaf <- t(mod_report$natage_f[,1:(spinup_years-1)])
                    naam <- t(mod_report$natage_m[,1:(spinup_years-1)])
                    naa_est[1:(spinup_years-1),,1,] <- naaf
                    naa_est[1:(spinup_years-1),,2,] <- naam

                    F_ll_f <- t(mod_report$F_ll_f[,1:(spinup_years-1)])
                    F_ll_m <- t(mod_report$F_ll_m[,1:(spinup_years-1)])
                    F_tw_f <- t(mod_report$F_trwl_f[,1:(spinup_years-1)])
                    F_tw_m <- t(mod_report$F_trwl_m[,1:(spinup_years-1)])

                    faa_est[1:(spinup_years-1),,1,,1] <- F_ll_f
                    faa_est[1:(spinup_years-1),,2,,1] <- F_ll_m
                    faa_est[1:(spinup_years-1),,1,,2] <- F_tw_f
                    faa_est[1:(spinup_years-1),,2,,2] <- F_tw_m

                }

                naa_est[y,,1,] <- mod_report$natage_f[,y]
                naa_est[y,,2,] <- mod_report$natage_m[,y]
                naa_proj <- naa_est[y,,,, drop=FALSE]

                faa_est[y,,1,,1] <- mod_report$F_ll_f[,y]
                faa_est[y,,2,,1] <- mod_report$F_ll_m[,y]
                faa_est[y,,1,,2] <- mod_report$F_trwl_f[,y]
                faa_est[y,,2,,2] <- mod_report$F_trwl_m[,y]
                prop_fs <- apply(faa_est[y,,,1,, drop=FALSE], 5, max)/sum(apply(faa_est[y,,,1,, drop=FALSE], 5, max))

                sel_est <- array(NA, dim=c(1, 30, 2, 1, 2))
                sel_est[1,,1,1,1] <- selex %>% filter(gear == "fixed", sex == "female", time_block == 3) %>% pull(value)
                sel_est[1,,1,1,2] <- selex %>% filter(gear == "trawl", sex == "female", time_block == 1) %>% pull(value)
                sel_est[1,,2,1,1] <- selex %>% filter(gear == "fixed", sex == "male", time_block == 3) %>% pull(value)
                sel_est[1,,2,1,2] <- selex %>% filter(gear == "trawl", sex == "male", time_block == 1) %>% pull(value)

                sel <- sel_est

                model_outs$mods[[(y+1)-spinup_years]] <- mod_out$model
                model_outs$fits[[(y+1)-spinup_years]] <- mod_out$opt
                model_outs$reps[[(y+1)-spinup_years]] <- mod_out$report

            }


            # Solve for reference points, F from the HCR,
            # and compute TAC for the next year. Note that
            # selectivity for RP calculations is weighted
            # by terminal year F.
            joint_selret <- calculate_joint_selret(sel, dp_y$ret, prop_fs)

            # reference points are all female based
            ref_pts <- calculate_ref_points(
                nages=nages,
                mort = dp_y$mort[,,1,],
                mat = dp_y$mat[,,1,],
                waa = dp_y$waa[,,1,],
                sel =  joint_selret$sel[,,1,,drop=FALSE],
                ret = joint_selret$ret[,,1,,drop=FALSE],
                avg_rec = mean(rec)/2,
                spr_target = mp$ref_points$spr_target
            )

            hcr_parameters <- list(ref_pts=ref_pts, naa=naa_proj, dem_params=dp_y)
            if(!is.na(mp$hcr$extra_pars)){
                hcr_parameters <- c(hcr_parameters, mp$hcr$extra_pars)
            }

            hcr_out <- do.call(mp$hcr$func, hcr_parameters)
            if(mp$hcr$units != "F"){
                hcr_out <- afscOM::find_F(
                    f_guess=0.10,
                    naa = naa_proj,
                    waa = dp_y$waa,
                    mort = dp_y$mort,
                    selex = joint_selret$sel,
                    ret = joint_selret$ret,
                    dmr = afscOM::subset_matrix(dp_y$dmr[,,,,1,drop=FALSE], 1, d=5, drop=TRUE),
                    prov_catch = hcr_out
                )
            }

            hcr_F[y] <- hcr_out
            #hcr_F[y] <- npfmc_tier3_F(assessment_ssb, ref_pts$B40, ref_pts$F40)

            mgmt_out <- simulate_TAC(
                hcr_F = hcr_F[y], 
                naa = naa_proj, 
                recruitment = mean(rec)/2, 
                joint_sel = joint_selret$sel, 
                dem_params = dp_y, 
                hist_tac = tac[y,1,1,1],
                hcr_options = mp$hcr$extra_options,
                options = mp$management
            )

            abc[y+1,1,1,1] <- mgmt_out$abc
            tac[y+1,1,1,1] <- mgmt_out$tac
            exp_land[y+1,1,1,1] <- mgmt_out$land
            landings[y+1] <- mgmt_out$land
            hcr_f[y,,,] <- hcr_F[y] 
        }

        if(y %% 10 == 0){
            print(paste0("Sim ", seed, ": ", y, "/", nyears_input, " complete."))
        }   
    }

    file.remove(paste0("data/sablefish_em_data_curr_",file_suffix,".RDS"))
    file.remove(paste0("data/sablefish_em_par_curr_",file_suffix,".RDS"))

    return(afscOM::listN(land_caa, disc_caa, caa, faa, faa_est, naa, naa_est, out_f, exp_land, hcr_f, abc, tac, survey_obs, model_outs))

}
