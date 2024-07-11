average_catch <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL){
    avg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
            as_tibble() %>%
            filter(time > 64) %>%
            group_by(time, sim, om, hcr) %>%
            mutate(total_catch = sum(value))
    
    if(!is.null(extra_filter)){
        avg_catch <- avg_catch %>% filter(eval(extra_filter))
    }
            
    return(
        avg_catch %>%
            group_by(om, hcr) %>%
            median_qi(total_catch, .width=interval_widths, .simple_names=FALSE)
    )
}

average_ssb <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL){
    
    agv_ssb <- get_ssb_biomass(model_runs, extra_columns) %>%
            filter(L1 != "naa_est", time > 64)

    if(!is.null(extra_filter)){
        agv_ssb <- agv_ssb %>% filter(eval(extra_filter))
    }
    
    return(
         agv_ssb %>%
            group_by(om, hcr) %>%
            median_qi(spbio, .width=interval_widths, .simple_names=FALSE)
    )
}

average_annual_catch_variation <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL){
    
    avg_catch_var <- get_landed_catch(model_runs, extra_columns) %>% 
            filter(time > 64)

    if(!is.null(extra_filter)){
        avg_catch_var <- avg_catch_var %>% filter(eval(extra_filter))
    }
    
    return(
         avg_catch_var %>%
            group_by(sim, hcr, om) %>%
            summarise(
                aav = aav(total_catch)
            ) %>%
            group_by(hcr, om) %>%
            median_qi(aav, .width=interval_widths, .simple_names=FALSE)
    )
}

average_proportion_catch_large <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL){
    
    prop_lg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
            as_tibble() %>%
            filter(time > 64) %>%
            mutate(
                size_group = case_when(
                    age < 5 ~ "Small",
                    age < 9 ~ "Medium",
                    TRUE ~ "Large"
                )
            ) %>%
            group_by(time, sim, size_group, om, hcr) %>%
            summarise(total_catch = sum(value)) %>%
            ungroup() %>%
            pivot_wider(names_from = "size_group", values_from="total_catch") %>%
            rowwise() %>%
            mutate(
                total_catch = sum(Large, Medium, Small)
            ) %>%
            mutate(across(Large:Small, ~ ./total_catch)) %>%
            select(-total_catch) %>%
            ungroup() %>%
            pivot_longer(Large:Small, names_to="size_group", values_to="catch")

    if(!is.null(extra_filter)){
        prop_lg_catch <- prop_lg_catch %>% filter(eval(extra_filter))
    }
    
    return(
         prop_lg_catch %>%
            group_by(size_group, om, hcr) %>%
            median_qi(catch, .width=interval_widths, .simple_names=FALSE)
    )
}

average_proportion_biomass_old <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL){
    
    prop_old_biomass <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
            as_tibble() %>%
            # join WAA and maturity-at-age for computing SSB
            left_join(
                melt(om_list$om1$dem_params$waa, value.name="weight"), 
                by=c("time", "age", "sex")
            ) %>%
            mutate(bio = value*weight) %>%
            select(time, age, sex, sim, om, hcr, bio) %>%
            filter(time > 64) %>%
            mutate(
                age_group = case_when(
                    age < 7 ~ "Young",
                    age < 21 ~ "Adult",
                    TRUE ~ "Old"
                )
            ) %>%
            group_by(time, sim, age_group, om, hcr) %>%
            summarise(total_bio = sum(bio)) %>%
            ungroup() %>%
            pivot_wider(names_from = "age_group", values_from="total_bio") %>%
            rowwise() %>%
            mutate(
                total_bio = sum(Young, Adult, Old)
            ) %>%
            mutate(across(Adult:Young, ~ ./total_bio)) %>%
            select(-total_bio) %>%
            ungroup() %>%
            pivot_longer(Adult:Young, names_to="age_group", values_to="bio") %>%
            group_by(age_group, om, hcr)

    if(!is.null(extra_filter)){
        prop_old_biomass <- prop_old_biomass %>% filter(eval(extra_filter))
    }
    
    return(
         prop_old_biomass %>%
            median_qi(bio, .width=interval_widths, .simple_names=FALSE) %>%
            arrange(.width, hcr, om, age_group)
    )
}

average_annual_value <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL){
    avg_rel_value <- bind_mse_outputs(model_runs, "caa", extra_columns2) %>%
        as_tibble() %>%
        filter(time > 64) %>%
        mutate(
            size_group = case_when(
                age < 5 ~ "Small",
                age < 9 ~ "Medium",
                TRUE ~ "Large"
            )
        ) %>%
        group_by(time, sim, size_group, om, hcr) %>%
        summarise(total_catch = sum(value)) %>%
        mutate(
            relative_value = case_when(
                size_group == "Large" ~ total_catch*1,
                size_group == "Medium" ~ total_catch*0.55619,
                size_group == "Small" ~ total_catch*0.33711
            )
        ) %>%
        group_by(time, sim, om, hcr) %>%
        summarise(total_value = sum(relative_value)) 

        if(!is.null(extra_filter)){
            avg_rel_value <- avg_rel_value %>% filter(eval(extra_filter))
        }
        
        return(
            avg_rel_value %>%
                group_by(sim, om, hcr) %>%
                summarise(annual_value = mean(total_value)) %>%
                group_by(om, hcr) %>%
                median_qi(annual_value, .width=interval_widths, .simple_names=FALSE)
        )
}

performance_metric_summary <- function(model_runs, extra_columns, interval_widths, extra_filter=NULL){
    # Average Catch Across Projection Period
    avg_catch <- average_catch(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter) %>% reformat_ggdist_long(n=2)

    # Average SSB Across Projection Period
    avg_ssb <- average_ssb(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter) %>% reformat_ggdist_long(n=2)

    # Average Annual Catch Variation Across Projection Period
    avg_variation <- average_annual_catch_variation(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter) %>% reformat_ggdist_long(n=2)

    # Average proportion of catch that is "large"
    avg_catch_lg <- average_proportion_catch_large(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter) %>% 
        filter(size_group == "Large") %>%
        select(-size_group) %>%
        reformat_ggdist_long(n=2)

    # Average proportion of population that is "old"
    avg_pop_old <- average_proportion_biomass_old(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter) %>% 
        filter(age_group == "Old") %>%
        select(-age_group) %>%
        reformat_ggdist_long(n=2)

    # Average annual value
    annual_value <- average_annual_value(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter) %>% reformat_ggdist_long(n=2)

    perf_data <- bind_rows(avg_catch, avg_ssb, avg_variation, avg_catch_lg, avg_pop_old, annual_value) %>%
        mutate(name=factor(
                        name, 
                        levels=c("total_catch", "spbio", "aav", "catch", "bio", "annual_value"), 
                        labels=c("Catch", "SSB", "Catch AAV", "Large Catch", "Old SSB", "Annual Value")
                    )
        )

    return(afscOM::listN(avg_catch, avg_ssb, avg_variation, avg_catch_lg, avg_pop_old, annual_value, perf_data))

}
