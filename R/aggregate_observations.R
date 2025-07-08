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

#' Generate Regionally Weighted Aggregated Age Compositions 
#' 
#' Create a single-region aggregate age composition from a spatially-explicit
#' OM. Regional age-compositions are weighted by modifying the ISS to the
#' random multinomial draw based on the proportion of stock biomass or total
#' landings in each region. This deviates slightly from the operational method,
#' where regional weights come from the LL RPW index and CAA in numbers of
#' inidviduals.
#' 
#' @param dem_params list of demographic parameters from the OM [1, nages, nsexes, nregions, nfleets]
#' @param naa array of numbers at age [1, nages, nsexes, nregions]
#' @param caa array of catch at age [1, nages, nsexes, nregions
#' @param obs_pars list of observation parameters from OM
#' 
#' @return array of age compositions [1, nages, nsexes, 1, nfleets+surveys]
#' 
#' @export generate_aggregate_comps
#' 
generate_aggregate_comps <- function(dem_params, naa, caa, obs_pars){
    nfleets <- dim(dem_params$sel)[5]
    nsurveys <- dim(dem_params$surv_sel)[5]

    tmp <- array(NA, dim=c(1, 30, 2, 1, nfleets+nsurveys))

    for(i in 1:(nfleets+nsurveys)){
        ISS <- obs_pars$ac_samps[i]
        agg_sex <- obs_pars$acs_agg_sex[i]
        as_int <- obs_pars$ac_as_integers[i]

        is_survey <- obs_pars$is_survey[i]
        if(is_survey){
            selex <- subset_matrix(dp_y$surv_sel, r=i-2, d=5, drop=TRUE)
            weights <- apply(naa*dp_y$waa, c(1, 4), sum)
        }else{
            selex <- subset_matrix(dp_y$sel, r=i, d=5, drop=TRUE)
            weights <- apply(caa, c(1, 4), sum) # this is landings (caa*weight)
        }

        agg_comp <- simulate_weighted_comps(
            ac = naa,
            weight_type = 1,
            weights = weights,
            selex = selex,
            total_samples = ISS,
            aggregate_sex = agg_sex,
            as_integers = as_int
        )

        tmp[,,,,i] <- agg_comp
    }

    return(tmp)
}