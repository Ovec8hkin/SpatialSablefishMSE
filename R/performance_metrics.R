average_catch <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
    avg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
            as_tibble() %>%
            filter(time > 64) %>%
            group_by(time, sim, om, hcr) %>%
            summarise(total_catch = sum(value))
    
    if(!is.null(relative)){
        avg_catch <- avg_catch %>%
            group_by(sim, om, hcr) %>%
            pivot_wider(names_from=hcr, values_from = total_catch) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(4:(ncol(.)), names_to="hcr", values_to="total_catch")
    }

    if(!is.null(extra_filter)){
        avg_catch <- avg_catch %>% filter(eval(extra_filter))
    }
            
    return(
        avg_catch %>%
            group_by(om, hcr) %>%
            median_qi(total_catch, .width=interval_widths, .simple_names=FALSE)
    )
}

total_catch <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=TRUE){
    tot_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
            as_tibble() %>%
            filter(time > 64) %>%
            group_by(sim, om, hcr) %>%
            mutate(total_catch = sum(value))
    
    if(!is.null(relative)){
        tot_catch <- tot_catch %>%
            group_by(sim, om, hcr) %>%
            pivot_wider(names_from=hcr, values_from = total_catch) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(4:(ncol(.)), names_to="hcr", values_to="total_catch")
    }

    if(!is.null(extra_filter)){
        tot_catch <- tot_catch %>% filter(eval(extra_filter))
    }
            
    return(
        tot_catch %>%
            group_by(om, hcr) %>%
            median_qi(total_catch, .width=interval_widths, .simple_names=FALSE)
    )
}

prop_years_catch <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), catch_threshold, extra_filter=NULL, relative=NULL){
    catch_years <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
            as_tibble() %>%
            filter(time > 64) %>%
            group_by(time, sim, om, hcr) %>%
            mutate(
                total_catch = sum(value),
            ) %>%
            group_by(sim, om, hcr) %>%
            summarise(
                num_years = sum(total_catch >= catch_threshold)
            )

    if(!is.null(relative)){
        catch_years <- catch_years %>%
            group_by(sim, om, hcr) %>%
            pivot_wider(names_from=hcr, values_from = num_years) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(4:(ncol(.)), names_to="hcr", values_to="num_years")
    }
    
    if(!is.null(extra_filter)){
        catch_years <- catch_years %>% filter(eval(extra_filter))
    }
            
    return(
        catch_years %>%
            group_by(om, hcr) %>%
            median_qi(num_years, .width=interval_widths, .simple_names=FALSE)
    )
}

average_ssb <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
    
    agv_ssb <- get_ssb_biomass(model_runs, extra_columns) %>%
            ungroup() %>%
            filter(L1 != "naa_est", time > 64) %>%
            select(-L1)

    if(!is.null(relative)){
        agv_ssb <- agv_ssb %>%
            group_by(sim, om, hcr) %>%
            pivot_wider(names_from=hcr, values_from = spbio) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(4:(ncol(.)), names_to="hcr", values_to="spbio")
    }

    if(!is.null(extra_filter)){
        agv_ssb <- agv_ssb %>% filter(eval(extra_filter))
    }
    
    return(
         agv_ssb %>%
            group_by(om, hcr) %>%
            median_qi(spbio, .width=interval_widths, .simple_names=FALSE)
    )
}

average_annual_catch_variation <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
    
    avg_catch_var <- get_landed_catch(model_runs, extra_columns) %>% 
            filter(time > 64) %>%
            group_by(sim, hcr, om) %>%
            summarise(
                aav = aav(total_catch)
            )

    if(!is.null(relative)){
        avg_catch_var <- avg_catch_var %>%
            group_by(sim, om, hcr) %>%
            pivot_wider(names_from=hcr, values_from = aav) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(4:(ncol(.)), names_to="hcr", values_to="aav")
    }

    if(!is.null(extra_filter)){
        avg_catch_var <- avg_catch_var %>% filter(eval(extra_filter))
    }
    
    return(
         avg_catch_var %>%
            group_by(hcr, om) %>%
            median_qi(aav, .width=interval_widths, .simple_names=FALSE)
    )
}

average_proportion_catch_large <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
    
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

    if(!is.null(relative)){
        prop_lg_catch <- prop_lg_catch %>%
            group_by(sim, om, hcr, size_group) %>%
            pivot_wider(names_from=hcr, values_from = catch) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(5:(ncol(.)), names_to="hcr", values_to="catch")
    }

    if(!is.null(extra_filter)){
        prop_lg_catch <- prop_lg_catch %>% filter(eval(extra_filter))
    }
    
    return(
         prop_lg_catch %>%
            group_by(size_group, om, hcr) %>%
            median_qi(catch, .width=interval_widths, .simple_names=FALSE)
    )
}

average_proportion_biomass_old <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
    
    prop_old_biomass <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
            as_tibble() %>%
            # join WAA and maturity-at-age for computing SSB
            left_join(
                melt(om_list[[1]]$dem_params$waa, value.name="weight"), 
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
            pivot_longer(Adult:Young, names_to="age_group", values_to="bio") 

    if(!is.null(relative)){
        prop_old_biomass <- prop_old_biomass %>%
            group_by(sim, om, hcr, age_group) %>%
            pivot_wider(names_from=hcr, values_from = bio) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(5:(ncol(.)), names_to="hcr", values_to="bio")
    }

    if(!is.null(extra_filter)){
        prop_old_biomass <- prop_old_biomass %>% filter(eval(extra_filter))
    }
    
    return(
         prop_old_biomass %>%
            group_by(age_group, om, hcr) %>%
            median_qi(bio, .width=interval_widths, .simple_names=FALSE) %>%
            arrange(.width, hcr, om, age_group)
    )
}

average_annual_value <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
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

        if(!is.null(relative)){
            avg_rel_value <- avg_rel_value %>%
                group_by(sim, om, hcr) %>%
                pivot_wider(names_from=hcr, values_from = total_value) %>%
                mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
                pivot_longer(4:(ncol(.)), names_to="hcr", values_to="total_value")
        }

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

average_annual_dynamic_value <- function(model_runs, extra_columns, interval_widths=c(0.50, 0.80), extra_filter=NULL, relative=NULL){
    
    compute_dynamic_value <- function(landings, min_price_age, max_price_age, breakpoints=c(15, 30)){
        if(landings < breakpoints[1]){
            return(max_price_age)
        }else if(landings >= breakpoints[1] & landings <= breakpoints[2]){
            return(
                min_price_age + (breakpoints[2]-landings)/(breakpoints[2]-breakpoints[1])*(max_price_age-min_price_age)
            )
        }else{
            return(min_price_age)
        }
    }

    price_age_f_low <- c(0.597895623, 1.320303448, 1.320303448, 1.856562267, 2.610111345, 2.610111345, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875)
    price_age_m_low <- c(0.597895623, 0.597895623, 1.320303448, 1.320303448, 1.856562267, 1.856562267, 1.856562267, 1.856562267, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345)
    price_data_low <- matrix(c(price_age_f_low, price_age_m_low), nrow=length(price_age_f_low), ncol=2)
    dimnames(price_data_low) <- list("age"=2:31, "sex"=c("F", "M"))

    price_age_f_max <- c(7.917460094, 8.40756497, 8.40756497, 9.944657109, 11.46480347, 11.46480347, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658)
    price_age_m_max <- c(7.917460094, 7.917460094, 8.40756497, 8.40756497, 9.944657109, 9.944657109, 9.944657109, 9.944657109, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347)
    price_data_max <- matrix(c(price_age_f_max, price_age_m_max), nrow=length(price_age_f_max), ncol=2)
    dimnames(price_data_max) <- list("age"=2:31, "sex"=c("F", "M"))

    dyn_value <- bind_mse_outputs(model_runs, c("land_caa"), extra_columns) %>%
        as_tibble() %>%
        group_by(time, sim, om, hcr) %>%
        mutate(tot_catch = sum(value)) %>%
        filter(time > 64, fleet == "Fixed") %>%
        left_join(
            reshape2::melt(price_data_low) %>% rename(min_price=value),
            by = c("age", "sex")
        ) %>%
        left_join(
            reshape2::melt(price_data_max) %>% rename(max_price=value),
            by = c("age", "sex")
        ) %>%
        rowwise() %>%
        mutate(
            dyn_price = compute_dynamic_value(tot_catch, min_price, max_price)
        ) %>%
        group_by(time, sim, om, hcr) %>%
        summarise(total_value = sum(dyn_price*value)) 
                     # only calculate value for fixed-gera fleet

    if(!is.null(relative)){
        dyn_value <- dyn_value %>%
            group_by(sim, om, hcr) %>%
            pivot_wider(names_from=hcr, values_from = total_value) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(relative)))) %>%
            pivot_longer(4:(ncol(.)), names_to="hcr", values_to="total_value")
    }

    if(!is.null(extra_filter)){
        avg_rel_value <- avg_rel_value %>% filter(eval(extra_filter))
    }

    return(
        dyn_value %>%
            group_by(sim, om, hcr) %>%
            summarise(dyn_annual_value = mean(total_value)) %>%
            group_by(om, hcr) %>%
            median_qi(dyn_annual_value, .width=interval_widths, .simple_names=FALSE)
    )

}

performance_metric_summary <- function(model_runs, extra_columns, dem_params, interval_widths, extra_filter=NULL, relative=NULL){
    # Average Catch Across Projection Period
    avg_catch <- average_catch(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter, relative=relative) %>% reformat_ggdist_long(n=2)

    # Average SSB Across Projection Period
    avg_ssb <- average_ssb(model_runs, extra_columns2, interval_widths, extra_filter=extra_filter, relative=relative) %>% reformat_ggdist_long(n=2)

    # Average Annual Catch Variation Across Projection Period
    avg_variation <- average_annual_catch_variation(model_runs, extra_columns2, interval_widths, relative=relative, extra_filter=extra_filter) %>% reformat_ggdist_long(n=2)

    # Average proportion of catch that is "large"
    avg_catch_lg <- average_proportion_catch_large(model_runs, extra_columns2, interval_widths, relative=relative, extra_filter=extra_filter) %>% 
        filter(size_group == "Large") %>%
        select(-size_group) %>%
        reformat_ggdist_long(n=2)

    # Average proportion of population that is "old"
    avg_pop_old <- average_proportion_biomass_old(model_runs, extra_columns2, interval_widths, relative=relative, extra_filter=extra_filter) %>% 
        filter(age_group == "Old") %>%
        select(-age_group) %>%
        reformat_ggdist_long(n=2)

    # Average annual value
    annual_value <- average_annual_value(model_runs, extra_columns2, interval_widths, relative=relative, extra_filter=extra_filter) %>% reformat_ggdist_long(n=2)

    perf_data <- bind_rows(avg_catch, avg_ssb, avg_variation, avg_catch_lg, avg_pop_old, annual_value) %>%
        mutate(name=factor(
                        name, 
                        levels=c("total_catch", "spbio", "aav", "catch", "bio", "annual_value"), 
                        labels=c("Catch", "SSB", "Catch AAV", "Large Catch", "Old SSB", "Annual Value")
                    )
        )

    return(afscOM::listN(avg_catch, avg_ssb, avg_variation, avg_catch_lg, avg_pop_old, annual_value, perf_data))

}
