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
