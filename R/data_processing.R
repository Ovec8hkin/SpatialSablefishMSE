get_ssb_biomass <- function(model_runs, extra_columns, dem_params){
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            # join WAA and maturity-at-age for computing SSB
            left_join(
                melt(sable_om$dem_params$waa, value.name="weight"), 
                by=c("time", "age", "sex")
            ) %>%
            left_join(
                melt(sable_om$dem_params$mat, value.name="maturity"), 
                by=c("time", "age", "sex")
            ) %>%
            drop_na() %>%
            # compute derived quantities
            mutate(
                biomass = value*weight,
                spbio = value*weight*maturity
            )
    )
}

get_fishing_mortalities <- function(model_runs, extra_columns){
    return(
        bind_mse_outputs(model_runs, c("faa", "faa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            group_by(time, fleet, sim, L1, hcr) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                F = max(value)
            ) %>%
            ungroup() %>%
            group_by(time, sim, L1, hcr) %>%
            # total F is the sum of fleet-based Fs
            mutate(
                total_F = sum(F)
            ) %>%
            ungroup()
    )
}

get_recruits <- function(model_runs, extra_columns){
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            filter(age == 2) %>%
            group_by(time, L1, hcr, sim) %>%
            summarise(rec=sum(value))
    )
}

#' Get Reference Points from MSE Simulations
#' 
#' Derive fishing mortality and biomass reference points
#' from completed MSE simulations.
#' 
#' Note that this function hasn't been tested when multiple
#' MSE simulations are present in the `model_runs` list.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param dem_params list of demographic parameter matrices
#' @param year the simualtion year to calculate reference 
#' points at
#'
#' @export get_reference_points
#'
#' @example
#' \dontrun{
#'      model_runs <- list(mse1)
#'      extra_columns <- list(hcr="hcr1")
#'      get_reference_points(model_runs, extra_columns, om$dem_params, nyears)
#' }
#'
get_reference_points <- function(model_runs, extra_columns, dem_params, year){

    avg_recruitment <- get_recruits(model_runs, extra_columns) %>%
        filter(time <= year) %>%
        group_by(time) %>%
        summarise(rec=median(rec)) %>%
        pull(rec) %>%
        mean

    prop_fs <- get_fishing_mortalities(model_runs, extra_columns) %>%
        filter(L1 != "faa_est", time <= year) %>%
        group_by(time, sim, fleet) %>%
        summarise(
            prop_f = F/total_F
        ) %>%
        pivot_wider(names_from = "fleet", values_from="prop_f") %>%
        ungroup() %>%
        summarise(across(Fixed:Trawl, mean)) %>%
        slice(1) %>%
        unlist(., use.names=FALSE)

    joint_sel_f <- apply(dem_params$sel[year,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum)/max(apply(dem_params$sel[year,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum))

    ref_pts <- calculate_ref_points(
        30,
        dem_params$mort[year,,1,1],
        dem_params$mat[year,,1,1],
        dem_params$waa[year,,1,1],
        joint_sel_f,
        dem_params$ret[year,,1,1,1],
        avg_recruitment/2,
        spr_target = 0.40
    )

    return(ref_pts)

}
