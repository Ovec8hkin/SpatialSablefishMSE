aggregate_indices <- function(survey_obs, y){

    out_dims1 <- dim(survey_obs$rpns)
    out_dims1[1] <- y
    out_dims1[4] <- 1

    aggregate_obs <- list(
        rpns = array(apply(survey_obs$rpns, c(1, 5), \(x) sum(x))[1:y,], dim=out_dims1),
        rpws = array(apply(survey_obs$rpws, c(1, 5), \(x) sum(x))[1:y,], dim=out_dims1)
    )

    return(aggregate_obs)

}


generate_aggregate_comps <- function(y, nfleets, weight_type=1){

    dp_y <- afscOM::subset_dem_params(om$dem_params, y, d=1, drop=FALSE)

    tmp <- array(NA, dim=c(1, 30, 2, 1, nfleets))
    for(f in 1:nfleets){
        ISS <- om$model_options$obs_pars$ac_samps[f]
        agg_sex <- om$model_options$obs_pars$acs_agg_sex[f]
        as_int <- om$model_options$obs_pars$ac_as_integers[f]

        is_survey <- om$model_options$obs_pars$is_survey[f]
        if(is_survey){
            selex <- subset_matrix(dp_y$surv_sel, r=f-2, d=5, drop=TRUE)
            weights <- apply(model_runs$naa[y,,,,drop=FALSE]*om$dem_params$waa[y,,,,drop=FALSE], 4, sum)
        }else{
            selex <- subset_matrix(dp_y$sel, r=f, d=5, drop=TRUE)
            weights <- apply(model_runs$caa[y,,,,,drop=FALSE], 4, sum)
        }

        test <- simulate_aggregated_comp(
            ac = model_runs$naa[y,,,,drop=FALSE],
            weight_type = weight_type,
            weights = weights,
            selex = selex,
            total_samples = ISS,
            aggregate_sex = agg_sex,
            as_integers = as_int
        )

        tmp[,,,,f] <- test
    } 
    return(tmp)
}