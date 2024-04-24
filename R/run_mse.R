#' Run Management Strategy Evaluation
#' #'
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
run_mse <- function(om, hcr, ..., nyears_input=NA, spinup_years=64, seed=1120, file_suffix=""){
   
    assessment <- dget("data/sablefish_assessment_2023.rdat")
   
    # Load OM parameters into global environment
    list2env(om, env=environment())

    # Load OM dimensions into global environment
    list2env(afscOM::get_model_dimensions(dem_params$sel), env=environment())

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

    survey_obs <- list(
        ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
        ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
        tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
        ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions)),
        fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions))
    )

    model_outs = list(
        mods = vector("list", length=c(nyears-spinup_years)),
        fits = vector("list", length=c(nyears-spinup_years)),
        reps = vector("list", length=c(nyears-spinup_years))
    )

    naa[1,,,] = init_naa

    set.seed(seed)
    recruitment <- assessment$natage.female[,1]*2
    projected_recruitment <- sample(recruitment, size=nyears-length(recruitment)+1, replace=TRUE)
    recruitment <- c(recruitment, projected_recruitment)

    for(y in 1:nyears){
        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        removals_input <- landings[y]
        fleet.props <- unlist(lapply(model_options$fleet_apportionment, \(x) x[y]))

        prev_naa <- naa[y,,,, drop = FALSE]
        out_vars <- project(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=prev_naa,
            recruitment=recruitment[y+1],
            fleet.props = fleet.props,
            options=model_options
        )

        # update state
        land_caa[y,,,,] <- out_vars$land_caa_tmp
        disc_caa[y,,,,] <- out_vars$disc_caa_tmp
        caa[y,,,,] <- out_vars$caa_tmp
        faa[y,,,,] <- out_vars$faa_tmp
        naa[y+1,,,] <- out_vars$naa_tmp
        out_f[y,,,] <- sum(out_vars$F_f_tmp[1,1,1,1,])
        hcr_f[y,,,] <- hcr_F[y]

        survey_obs$ll_rpn[y,,,] <- out_vars$surv_obs$ll_rpn
        survey_obs$ll_rpw[y,,,] <- out_vars$surv_obs$ll_rpw
        survey_obs$tw_rpw[y,,,] <- out_vars$surv_obs$tw_rpw
        survey_obs$ll_acs[y,,,] <- out_vars$surv_obs$ll_ac_obs
        survey_obs$fxfish_acs[y,,,] <- out_vars$surv_obs$fxfish_caa_obs

        
        if((y+1) > spinup_years){

            naa_proj <- out_vars$naa_tmp
            rec <- recruitment[1:y]
            sel <- dp.y$sel
            prop_fs <- apply(out_vars$faa_tmp[,,,1,, drop=FALSE], 5, max)/sum(apply(out_vars$faa_tmp[,,,1,, drop=FALSE], 5, max))

            if(model_options$run_estimation){
                # Do all of the data formatting and running
                # of the TMB Sablefish model
                assess_inputs <- format_em_data(
                    nyears = y,
                    dem_params = dp.y,
                    land_caa = out_vars$land_caa_tmp,
                    survey_indices = out_vars$surv_obs,
                    fxfish_caa_obs = survey_obs$fxfish_acs[y-1,,,,drop=FALSE], # Age comp data is one year delayed
                    ll_ac_obs = survey_obs$ll_acs[y-1,,,,drop=FALSE], # Age comp data is one year delayed
                    model_options = model_options,
                    added_years = 1,
                    file_suffix = file_suffix
                )

                mod_out <- fit_TMB_model(assess_inputs$new_data, assess_inputs$new_parameters)  
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
            joint_self <- apply(sel[,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum)/max(apply(sel[,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum))
            joint_selm <- apply(sel[,,2,,,drop=FALSE]*prop_fs, c(1, 2), sum)/max(apply(sel[,,2,,,drop=FALSE]*prop_fs, c(1, 2), sum))
            joint_ret <- apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum))
            
            # reference points are all female based
            ref_pts <- calculate_ref_points(
                nages=nages,
                mort = dp.y$mort[,,1,],
                mat = dp.y$mat[,,1,],
                waa = dp.y$waa[,,1,],
                sel =  joint_self,
                ret = joint_ret,
                avg_rec = mean(rec)/2
            )

            hcr_F[y] <- match.fun(hcr)(ref_pts, naa_proj, dp.y, ...)
            #hcr_F[y] <- npfmc_tier3_F(assessment_ssb, ref_pts$B40, ref_pts$F40)

            joint_sel <- array(NA, dim=dim(out_vars$naa_tmp))
            joint_sel[,,1,] <- joint_self
            joint_sel[,,2,] <- joint_selm

            mgmt_out <- simulate_TAC(hcr_F[y], naa_proj, mean(rec)/2, joint_sel, dp.y, tac_land_reduction = 0.85)
            abc[y+1,1,1,1] <- mgmt_out$abc
            tac[y+1,1,1,1] <- mgmt_out$tac
            exp_land[y+1,1,1,1] <- mgmt_out$land
            landings[y+1] <- mgmt_out$land 
        }

        if(y %% 10 == 0){
            print(paste0("Sim ", seed, ": ", y, "/", nyears_input, " complete."))
        }   
    }

    file.remove(paste0("data/sablefish_em_data_curr_",file_suffix,".RDS"))
    file.remove(paste0("data/sablefish_em_par_curr_",file_suffix,".RDS"))

    return(afscOM::listN(land_caa, disc_caa, caa, faa, faa_est, naa, naa_est, out_f, exp_land, hcr_f, abc, tac, survey_obs, model_outs))

}
