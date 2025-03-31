#' Run Management Strategy Evaluation
#' 
#' Run an MSE simulation loop.
#'
#' @param om an operating model list object (should minimally contain a
#' list object with demographic parameters, an array of initial NAAs, a
#' list object specifying recruitment, and a list of model_options).
#' @param mp a list object specifying the management procedure to be applied during
#' the projection period
#' @param ... parameters to pass to the `hcr` function
#' @param run_estimation whether to run the estimation procedure or not. If FALSE, 
#' management procedures will use OM outputs rather than EM estimates to compute
#' future harvest from the management procedure function.
#' @param nyears_input number of years to simulate forward (will override the dimensions
#' of the demagraphic parameters matrices defined in the OM)
#' @param spinup_years number of years before estimation process should begin
#' @param seed random seed
#'
#' @export run_mse
#'
#' @example
#'
run_mse <- function(om, mp, mse_options, nyears_input=NA, seed=1120, file_suffix=""){
   
    spinup_years <- mse_options$n_spinup_years
    nyears_input <- mse_options$n_proj_years + mse_options$n_spinup_years

    # Setup what years to perform assessment in based on assessment_frequency
    # input. If input as a vector, use the vector literally. If input as a single
    # number, assumes an annual frequency.
    do_assessment <- generate_annual_frequency(mp$assessment_frequency, nyears_input)
    do_assessment[spinup_years] <- 1
    assessment_years <- which(do_assessment == 1)
    assessment_years[length(assessment_years)] <- nyears_input+1

    # Setup what years to perform surveys in based on *_survey_frequency
    # inputs. If input as a vector, use the vector literally. If input as a single
    # number, assumes an annual frequency.
    do_survey_ll <- generate_annual_frequency(mp$ll_survey_frequency, nyears_input - spinup_years)
    do_survey_tw <- generate_annual_frequency(mp$tw_survey_frequency, nyears_input - spinup_years)

    spatial_assessment <- dget(file.path(here::here(), "data", "spatial_sablefsh_inputs.rdat"))
   
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
        "region" = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
        "fleet" = c("Fixed", "Trawl")
    )

    landings <- array(0, dim=c(nyears+1, nfleets, nregions), dimnames=list("time"=1:(nyears+1), "fleet"=c("Fixed", "Trawl"), "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA")))
    hcr_F <- rep(0, nyears)
    out_f <- rep(0, nyears) # vector to store outputted F
    landings[1:spinup_years,,] <- aperm(spatial_assessment$catch[,1:spinup_years,], c(2, 3, 1))


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
    faa_est     = array(NA, dim=c(nyears, nages, nsexes, 1, nfleets), dimnames=list("time"=1:(nyears),   "age"=2:31, "sex"=c("F", "M"), "region"=c("alaska"), "fleet"=c("Fixed", "Trawl")))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA")))
    naa_est     = array(NA, dim=c(nyears,   nages, nsexes, 1), dimnames=list("time"=1:(nyears),   "age"=2:31, "sex"=c("F", "M"), "region"=c("alaska")))

    abc         = array(NA, dim=c(nyears+1, 1, 1, 1), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska"))
    tac         = array(NA, dim=c(nyears+1, 1, 1, 1), dimnames=list("time"=1:(nyears+1), 1, 1, "region"="Alaska"))
    exp_land    = array(NA, dim=c(nyears+1, 1, 1, nregions), dimnames=list("time"=1:(nyears+1), 1, 1, "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA")))
    hcr_f       = array(NA, dim=c(nyears+1,   1, 1, 1), dimnames=list("time"=1:(nyears+1),     1, 1, "region"="Alaska"))
    # out_f       = array(NA, dim=c(nyears,   1, 1, 1), dimnames=list("time"=1:nyears,     1, 1, "region"="Alaska"))

    f           = array(NA, dim=c(nyears, 1, 1, nregions, nfleets))
    global_rec_devs    = array(NA, dim=c(mse_options$n_proj_years, 1, 1, nregions), dimnames=list("time"=1:(mse_options$n_proj_years), 1, 1, "region"="Alaska"))

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
    hist_recruitment <- hist_recruitment[1:mse_options$recruitment_start]
    projected_recruitment <- do.call(recruitment$func, c(recruitment$pars, list(seed=seed)))
    if(!is.function(projected_recruitment)){
        full_recruitment <- c(hist_recruitment, projected_recruitment)
    }else{
        full_recruitment <- rep(NA, nyears)
        full_recruitment[1:length(hist_recruitment)] <- hist_recruitment
        set.seed(seed)
        rec_devs <- rlnorm(mse_options$n_proj_years, meanlog = 0, sdlog = 1.20)
        global_rec_devs[1:mse_options$n_proj_years, 1, 1, 1] <- rec_devs
    }
    
    for(y in 1:nyears_input){
        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp_y <- afscOM::subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        removals_input <- afscOM::subset_matrix(landings, y, d=1, drop=FALSE)
        # fleet.props <- subset_matrix(model_options$fleet_apportionment, y, d=1, drop=FALSE)
        # removals_input <- array(removals_input*fleet.props, dim=c(1, nfleets, nregions))
        # # region_props <- as.matrix(1)
        # rec_props <- as.matrix(1)

        # will work for any recruitment function that only requires
        # ssb as a yearly input (beverton holt and ricker should work fine)
        if(is.na(full_recruitment[y+1])){
            ssb <- sum(naa[y,,1,,drop=FALSE]*dp_y$waa[,,1,]*dp_y$mat[,,1,])
            full_recruitment[y+1] <- projected_recruitment(ssb, y-spinup_years+1) 
            full_recruitment[y+1] <- full_recruitment[y+1]*rec_devs[y-spinup_years+1] 
        }

        prev_naa <- naa[y,,,, drop = FALSE]
        out_vars <- project_single(
            removals = removals_input,
            dem_params=dp_y,
            prev_naa=prev_naa,
            recruitment=full_recruitment[y+1,],
            # fleet_props = fleet.props,
            # region_props = region_props,
            # rec_props = rec_props,
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

        
        if((y+1) > spinup_years && do_assessment[y]){

            naa_proj <- out_vars$naa_tmp
            rec <- full_recruitment[1:y]
            sel <- dp_y$sel
            prop_fs <- apply(out_vars$faa_tmp[,,,1,, drop=FALSE], 5, max)/sum(apply(out_vars$faa_tmp[,,,1,, drop=FALSE], 5, max))
            bio <- sum(naa_proj*dp_y$waa, na.rm=TRUE)

            if(mse_options$run_estimation && bio > 1){
                # Do all of the data formatting and running
                # of the TMB Sablefish model

                # First aggregate the observation data for use in a single-area assessment
                aggregated_survey_obs <- aggregate_observations(survey_obs, y)

                aggregate_land_caa <- array(apply(land_caa, c(1, 2, 3, 5), sum), dim=c(nyears, nages, nsexes, 1, nfleets))

                assess_inputs <- simulate_em_data_sex_disaggregate(
                    nyears = y,
                    dem_params = afscOM::subset_dem_params(om$dem_params, 1:y, d=1, drop=FALSE),
                    land_caa = aggregate_land_caa[1:y,,,,,drop=FALSE],
                    survey_indices = afscOM::subset_dem_params(aggregated_survey_obs, 1:y, d=1, drop=FALSE),
                    fxfish_caa_obs = afscOM::subset_matrix(aggregated_survey_obs$acs[1:y,,,,1,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    twfish_caa_obs = afscOM::subset_matrix(aggregated_survey_obs$acs[1:y,,,,2,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    ll_ac_obs = afscOM::subset_matrix(aggregated_survey_obs$acs[1:y,,,,3,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    tw_ac_obs = afscOM::subset_matrix(aggregated_survey_obs$acs[1:y,,,,4,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                    ll_srv_indic = do_survey_ll[(spinup_years:y)-spinup_years+1],
                    tw_srv_indic = do_survey_tw[(spinup_years:y)-spinup_years+1],
                    model_options = om$model_options,
                    added_years = y-spinup_years+1,
                    file_suffix = y
                )


                mod_out <- fit_TMB_model(
                    data = assess_inputs$new_data, 
                    parameters = assess_inputs$new_parameters,
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

            # NOTE: All demographic matrices are substted to the first spatial regions
            # by default. This does not matter so long as there are not regionally varying
            # demographic rates. If regionally variable demographics, needs to come up
            # with better way to do reference point calculations and projections.
            n_proj_years <- min(assessment_years[assessment_years > y]) - y - 1
            for(y2 in y:(y+n_proj_years)){
                # Solve for reference points, F from the HCR,
                # and compute TAC for the next year. Note that
                # selectivity for RP calculations is weighted
                # by terminal year F.
                joint_selret <- calculate_joint_selret(sel, dp_y$ret, prop_fs)

                # reference points are all female based
                ref_pts <- calculate_ref_points(
                    nages=nages,
                    mort = dp_y$mort[,,1,1],
                    mat = dp_y$mat[,,1,1],
                    waa = dp_y$waa[,,1,1],
                    sel =  joint_selret$sel[,,1,1,drop=FALSE],
                    ret = joint_selret$ret[,,1,1,drop=FALSE],
                    avg_rec = mean(rec)/2,
                    spr_target = mp$ref_points$spr_target
                )

                naa_proj[is.na(naa_proj)] <- 0
                hcr_parameters <- list(ref_pts=ref_pts, naa=naa_proj, dem_params=dp_y, avgrec=mean(rec))
                if(all(!is.na(mp$hcr$extra_pars))){
                    hcr_parameters <- c(hcr_parameters, mp$hcr$extra_pars)
                }

                hcr_out <- do.call(mp$hcr$func, hcr_parameters)
                if(mp$hcr$units != "F"){
                    hcr_out <- afscOM::find_F(
                        f_guess=0.10,
                        naa = naa_proj,
                        waa = afscOM::subset_matrix(dp_y$waa, 1, d=4, drop=FALSE),
                        mort = afscOM::subset_matrix(dp_y$mort, 1, d=4, drop=FALSE),
                        selex = joint_selret$sel,
                        ret = joint_selret$ret[,,,1, drop=FALSE],
                        dmr = afscOM::subset_matrix(afscOM::subset_matrix(dp_y$dmr[,,,,1,drop=FALSE], 1, d=5, drop=TRUE), 1, d=4, drop=FALSE),
                        prov_catch = hcr_out
                    )
                }

                # hcr_F[y+1] <- hcr_out
                #hcr_F[y] <- npfmc_tier3_F(assessment_ssb, ref_pts$B40, ref_pts$F40)
                
                mgmt_out <- simulate_TAC(
                    hcr_F = hcr_out, 
                    naa = naa_proj, 
                    recruitment = mean(rec)/2, 
                    joint_sel = joint_selret$sel, 
                    dem_params = afscOM::subset_dem_params(dp_y, 1, d=4, drop=FALSE),
                    hist_abc = abc[y2,1,1,1],
                    hcr_options = mp$hcr$extra_options,
                    options = mp$management
                )

                abc[y2+1,1,1,1] <- mgmt_out$abc
                tac[y2+1,1,1,1] <- mgmt_out$tac
                exp_land[y2+1,1,1,1] <- mgmt_out$land
                # This is default apportionment for now
                # TODO: generalize catch apportionment
                landings[y2+1,,] <- mgmt_out$land/nfleets/nregions
                hcr_f[y2+1,,,] <- hcr_out

                naa_proj <- mgmt_out$proj_N_new
            }
   
        }

        if(y %% 10 == 0){
            print(paste0("OM: ", om$name, "; MP: ", mp$name, "; Sim ", seed, ": ", y, "/", nyears_input, " complete."))
        }   
    }

    file.remove(paste0("data/sablefish_em_data_curr_",file_suffix,".RDS"))
    file.remove(paste0("data/sablefish_em_par_curr_",file_suffix,".RDS"))

    return(afscOM::listN(land_caa, disc_caa, caa, faa, faa_est, naa, naa_est, out_f, exp_land, hcr_f, abc, tac, global_rec_devs, survey_obs, model_outs))

}

