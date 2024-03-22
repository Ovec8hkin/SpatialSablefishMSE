#' Run Management Strategy Evaluation
#' #'
#' Run an MSE simulation loop.
#'
#' @param om path to a saved OM .RDS file (should contain
#' demagraphic parameters list and model options list)
#' @param hcr a function to compute the allowable F in the next year
#' @param ... parameters to pass to the `hcr` function
#' @param nsims number of random simulation to run
#' @param seed random seed
#'
#' @export run_mse
#'
#' @example
#'
run_mse <- function(om, hcr, ..., nsims=10, seed=1120){
   
    assessment <- dget("data/sablefish_assessment_2023.rdat")
   
    # Load OM parameters into global environment
    list2env(om, globalenv())

    # Load OM dimensions into global environment
    list2env(afscOM::get_model_dimensions(dem_params$sel), globalenv())

    dimension_names <- list(
        "time" = 1:nyears,
        "age"  = 2:31,
        "sex"  = c("F", "M"),
        "region" = "alaska",
        "fleet" = c("Fixed", "Trawl"),
        "sims" = 1:nsims
    )

    TACs <- rep(0, nyears)
    hcr_F <- rep(0, nyears)
    out_f <- rep(0, nyears) # vector to store outputted F
    TACs[1:64] <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])

    #' 6. Setup empty array to collect derived quantities from the OM
    #' It is left to the user to decide what information to store, and
    #' in what format they would like to store it.
    #'
    #' Here, we are storing all OM outputs as they are returned from the
    #' OM.
    land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
    tac         = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sims"=1:nsims))
    hcr_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sims"=1:nsims))
    out_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sims"=1:nsims))
    naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=1:nsims))
    naa[1,,,,] = init_naa

    set.seed(seed)
    seeds <- sample(1:(nsims*1000), size=nsims, replace=FALSE)
    # ref_as <- readRDS("data/agestruct_f40.RDS")
    # ref_naa <- ref_as$naa
    for(s in 1:nsims){
        print(s)
        set.seed(seeds[s])
        recruitment <- assessment$natage.female[,1]*2
        projected_recruitment <- sample(recruitment, size=nyears-length(recruitment)+1, replace=TRUE)
        recruitment <- c(recruitment, projected_recruitment)

        for(y in 1:(nyears-1)){

            # Subset the demographic parameters list to only the current year
            # and DO NOT drop lost dimensions.
            dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
            removals_input <- TACs[y]
            fleet.props <- unlist(lapply(model_options$fleet_apportionment, \(x) x[y]))

            prev_naa <- afscOM::subset_matrix(naa[y,,,,s, drop = FALSE], 1, 5, drop=TRUE)
            out_vars <- project(
                removals = removals_input,
                dem_params=dp.y,
                prev_naa=prev_naa,
                recruitment=recruitment[y+1],
                fleet.props = fleet.props,
                options=model_options
            )

            # update state
            land_caa[y,,,,,s] <- out_vars$land_caa_tmp
            disc_caa[y,,,,,s] <- out_vars$disc_caa_tmp
            caa[y,,,,,s] <- out_vars$caa_tmp
            faa[y,,,,,s] <- out_vars$faa_tmp
            naa[y+1,,,,s] <- out_vars$naa_tmp
            out_f[y,,,,s] <- sum(out_vars$F_f_tmp[1,1,1,1,])
            tac[y,,,,s] <- TACs[y]
            hcr_f[y,,,,s] <- hcr_F[y]

            #new_F <- 0.085
            if((y+1) > 64){
                joint_self <- apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum))
                joint_selm <- apply(dp.y$sel[,,2,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum))
                joint_ret <- apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum))
                
                # reference points are all female based
                ref_pts <- calculate_ref_points(
                    nages=nages,
                    mort = dp.y$mort[,,1,],
                    mat = dp.y$mat[,,1,],
                    waa = dp.y$waa[,,1,],
                    sel =  joint_self,
                    ret = joint_ret,
                    avg_rec = mean(recruitment)/2
                )

                ssb <- apply(out_vars$naa_tmp[,,1,]*dp.y$waa[,,1,,drop=FALSE]*dp.y$mat[,,1,,drop=FALSE], 1, sum)
                # hcr_F[y] <- as_scalar_threshold_f(
                #                 ssb/ref_pts$B40, 
                #                 naa=out_vars$naa_tmp[,,1,], 
                #                 ref_naa=ref_naa,
                #                 as_func = shannon_diversity,
                #                 #ages = 2:31,
                #                 f_min=0,
                #                 f_max=ref_pts$F40,
                #                 lrp=0.05,
                #                 urp=1.0
                #             )

                hcr_F[y] <- match.fun(hcr)(ref_pts, out_vars$naa_tmp, dp.y, ...)

                joint_sel <- array(NA, dim=dim(out_vars$naa_tmp))
                joint_sel[,,1,] <- joint_self
                joint_sel[,,2,] <- joint_selm

                TACs[y+1] <- simulate_TAC(hcr_F[y], out_vars$naa_tmp, mean(recruitment)/2, joint_sel, dp.y)
            }   
        }
    }

    return(afscOM::listN(naa, caa, faa, out_f))

}
