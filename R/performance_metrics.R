#' Compute Average Catch across projection period
#' 
#' Compute average catch (median and CIs) per year across all years and simulation seeds, 
#' for each combination of operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_catch
#'
#' @example
#'
average_catch <- function(
    model_runs, 
    extra_columns, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA),  
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){

    group_columns <- c("sim", summarise_by)
    
    avg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_times(time_horizon) %>%
        group_by(across(all_of(c("time", group_columns)))) %>%
        summarise(annual_catch = sum(value)) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "annual_catch",
            rel_value = relative,
            grouping=group_columns
        )
    
    if(!is.null(extra_filter)){
        avg_catch <- avg_catch %>% filter(eval(extra_filter))
    }
            
    return(
        avg_catch %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(annual_catch, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute Toatl Catch over projection period
#' 
#' Compute total catch (median and CIs) over all years and simulation seeds, 
#' for each combination of operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export total_catch
#'
#' @example
#'
total_catch <- function(
    model_runs, 
    extra_columns, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)
    tot_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_times(time_horizon = time_horizon) %>%
        group_by(across(all_of(group_columns))) %>%
        summarise(total_catch = sum(value)) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "total_catch",
            rel_value = relative,
            grouping=group_columns
        )

    if(!is.null(extra_filter)){
        tot_catch <- tot_catch %>% filter(eval(extra_filter))
    }
            
    return(
        tot_catch %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(total_catch, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute proportion of years catch exceeds a threshold level across projection period
#' 
#' Compute average proportion of years where catch exceeds a threshold level
#' (median and CIs) per year across all years and simulation seeds, 
#' for each combination of operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export prop_years_catch
#'
#' @example
#'
prop_years_catch <- function(
    model_runs, 
    extra_columns, 
    catch_threshold, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)
    catch_years <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_times(time_horizon = time_horizon) %>%
        group_by(across(all_of(c("time", group_columns)))) %>%
        summarise(
            total_catch = sum(value),
        ) %>%
        group_by(across(all_of(group_columns))) %>%
        summarise(
            num_years = sum(total_catch >= catch_threshold)/n()
        ) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "num_years",
            rel_value = relative,
            grouping = group_columns
        )

    if(!is.null(extra_filter)){
        catch_years <- catch_years %>% filter(eval(extra_filter))
    }
            
    return(
        catch_years %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(num_years, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute Average SSB across projection period
#' 
#' Compute average SSB (median and CIs) per year across all years and 
#' simulation seeds, for each combination of operating models and management 
#' procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param dem_params demographic parameters matrices from OM
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_ssb
#'
#' @example
#'
average_ssb <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)

    agv_ssb <- get_ssb_biomass(model_runs, extra_columns, dem_params) %>%
            ungroup() %>%
            filter(L1 != "naa_est") %>%
            select(-L1) %>%
            filter_times(time_horizon=time_horizon) %>%
            relativize_performance(
                rel_column = "hcr",
                value_column = "spbio",
                rel_value = relative,
                grouping = group_columns
            )
            
    if(!is.null(extra_filter)){
        agv_ssb <- agv_ssb %>% filter(eval(extra_filter))
    }

    return(
         agv_ssb %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(spbio, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute Average Age across projection period
#' 
#' Compute average age (median and CIs) per year across all years and 
#' simulation seeds, for each combination of operating models and management 
#' procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param dem_params demographic parameters matrices from OM
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_age
#'
#' @example
#'
average_age <- function(
    model_runs, 
    extra_columns,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)

    avg_age <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
            as_tibble() %>%
            ungroup() %>%
            group_by(time, age, sim, om, hcr) %>%
            mutate(value = sum(value)) %>%
            filter(sex == "F") %>%
            ungroup() %>%
            group_by(time, sim, hcr, om) %>%
            summarise(
                avg_age = compute_average_age(value, 2:31)
            ) %>%
            filter_times(time_horizon=time_horizon) %>%
            relativize_performance(
                rel_column = "hcr",
                value_column = "avg_age",
                rel_value = relative,
                grouping = group_columns
            )
            
    if(!is.null(extra_filter)){
        avg_age <- avg_age %>% filter(eval(extra_filter))
    }

    return(
         avg_age %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(avg_age, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute Average ABI across projection period
#' 
#' Compute average ABI (median and CIs; Griffiths et al. 2024) per year 
#' across all years and simulation seeds, for each combination of operating 
#' models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param ref_naa reference age structure for ABI computation
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_abi
#'
#' @example
#'
average_abi <- function(
    model_runs, 
    extra_columns,
    ref_naa,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)

    avg_abi <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
            as_tibble() %>%
            ungroup() %>%
            group_by(time, age, sim, om, hcr) %>%
            mutate(value = sum(value)) %>%
            filter(sex == "F") %>%
            ungroup() %>%
            group_by(time, sim, hcr, om) %>%
            summarise(
                avg_abi = abi(value, ref_naa)
            ) %>%
            filter_times(time_horizon=time_horizon) %>%
            relativize_performance(
                rel_column = "hcr",
                value_column = "avg_abi",
                rel_value = relative,
                grouping = group_columns
            )
            
    if(!is.null(extra_filter)){
        avg_abi <- avg_abi %>% filter(eval(extra_filter))
    }

    return(
         avg_abi %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(avg_abi, .width=interval_widths, .simple_names=FALSE)
    )
}



#' Compute Average Annual Catch Variation (AAV) across projection period
#' 
#' Compute average annual catch variation (median and CIs) per year 
#' across all years and simulation seeds, for each combination of 
#' operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_annual_catch_variation
#'
#' @example
#'
average_annual_catch_variation <- function(
    model_runs, 
    extra_columns,
    interval_widths=c(0.50, 0.80),
    time_horizon = c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)
    avg_catch_var <- get_landed_catch(model_runs, extra_columns) %>%
        filter_times(time_horizon = time_horizon) %>%
        group_by(across(all_of(group_columns))) %>%
        summarise(
            aav = aav(total_catch)
        ) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "aav",
            rel_value = relative,
            grouping = group_columns
        )

    if(!is.null(extra_filter)){
        avg_catch_var <- avg_catch_var %>% filter(eval(extra_filter))
    }

    return(
         avg_catch_var %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(aav, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute average proportion of catch that is "Large" across projection period
#' 
#' Compute average proportion of catch that is "large" (median and CIs) 
#' per year across all years and simulation seeds, for each combination 
#' of operating models and management procedures. "Large" fish are those 
#' >9yo.
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_proportion_catch_large
#'
#' @example
#'
average_proportion_catch_large <- function(
    model_runs, 
    extra_columns, 
    interval_widths=c(0.50, 0.80),
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)
    prop_lg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_times(time_horizon = time_horizon) %>%
        mutate(
            size_group = case_when(
                age < 5 ~ "Small",
                age < 9 ~ "Medium",
                TRUE ~ "Large"
            )
        ) %>%
        group_by(across(all_of(c("time", "size_group", group_columns)))) %>%
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
        pivot_longer(Large:Small, names_to="size_group", values_to="catch") %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "catch",
            rel_value = relative,
            grouping = c("size_group", group_columns)
        )

    if(!is.null(extra_filter)){
        prop_lg_catch <- prop_lg_catch %>% filter(eval(extra_filter))
    }
    
    return(
         prop_lg_catch %>%
            filter(size_group == "Large") %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(catch, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute average proportion of population that is "Old" across projection period
#' 
#' Compute aaverage proportion of population that is "Old" (median and CIs) 
#' per year across all years and simulation seeds, for each combination 
#' of operating models and management procedures. "Old" fish are those 
#' >21yo.
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_proportion_biomass_old
#'
#' @example
#'
average_proportion_biomass_old <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
    group_columns <- c("sim", summarise_by)
    prop_old_biomass <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
        as_tibble() %>%
        # join WAA and maturity-at-age for computing SSB
        left_join(
            melt(dem_params$waa, value.name="weight"), 
            by=c("time", "age", "sex")
        ) %>%
        mutate(bio = value*weight) %>%
        filter_times(time_horizon = time_horizon) %>%
        mutate(
            age_group = case_when(
                age < 7 ~ "Young",
                age < 21 ~ "Adult",
                TRUE ~ "Old"
            )
        ) %>%
        group_by(across(all_of(c("time", "age_group", group_columns)))) %>%
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
        relativize_performance(
            rel_column = "hcr",
            value_column = "bio",
            rel_value = relative,
            grouping = c("age_group", group_columns)
        )

    if(!is.null(extra_filter)){
        prop_old_biomass <- prop_old_biomass %>% filter(eval(extra_filter))
    }
    
    return(
         prop_old_biomass %>%
            filter(age_group == "Old") %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(bio, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute Average Annual Value of Catch across projection period
#' 
#' Compute average annual value of landed catch (median and CIs) 
#' per year across all years and simulation seeds, for each combination 
#' of operating models and management procedures. Fish >9yo are considered
#' "high value" and are assigned a relative value of 1, fish 5-9yo are
#' considered "medium value" and are asisgned a relative value of 0.556, 
#' and fish <5yo are considered "low value" and are assigned a relative
#' valye of 0.337.
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_annual_value
#'
#' @example
#'
average_annual_value <- function(
    model_runs, 
    extra_columns, 
    interval_widths=c(0.50, 0.80),
    time_horizon = c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){

    group_columns <- c("sim", summarise_by)
    avg_rel_value <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_times(time_horizon = time_horizon) %>%
        mutate(
            size_group = case_when(
                age < 5 ~ "Small",
                age < 9 ~ "Medium",
                TRUE ~ "Large"
            )
        ) %>%
        group_by(across(all_of(c("time", "size_group", group_columns)))) %>%
        summarise(total_catch = sum(value)) %>%
        mutate(
            relative_value = case_when(
                size_group == "Large" ~ total_catch*1,
                size_group == "Medium" ~ total_catch*0.55619,
                size_group == "Small" ~ total_catch*0.33711
            )
        ) %>%
        group_by(across(all_of(c("time", group_columns)))) %>%
        summarise(total_value = sum(relative_value)) %>%
        group_by(across(all_of(group_columns))) %>%
        summarise(annual_value = mean(total_value)) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "annual_value",
            rel_value = relative,
            grouping = group_columns
        )

    if(!is.null(extra_filter)){
        avg_rel_value <- avg_rel_value %>% filter(eval(extra_filter))
    }

    return(
        avg_rel_value %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(annual_value, .width=interval_widths, .simple_names=FALSE)
    )
}

#' Compute Average Annual Dynamic Value of Catch across projection period
#' 
#' Compute average annual dynamic value of landed catch from the
#' fixed-gear fleet (median and CIs) per year across all years and 
#' simulation seeds, for each combination of operating models and management 
#' procedures. Value is computed based on an assumed linear relationship
#' between catch and price (assuming that increased catch floods existing
#' markets and decreases prices due to market abundance), where prices are 
#' constant at a "maximum" level when catch <15k mt, and are constant at a 
#' "minimum" level when catch >30k mt. Pricing data from June 2024 NPFMC
#' Small Sablefish Release Analyses.
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export average_annual_dynamic_value
#'
#' @example
#'
average_annual_dynamic_value <- function(
    model_runs, 
    extra_columns, 
    interval_widths=c(0.50, 0.80),
    time_horizon = c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    
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
    
    max_price <- max(c(price_age_f_max, price_age_m_max))
    price_data_low <- price_data_low/price_data_max
    price_data_max <- price_data_max/max_price
    price_data_low <- price_data_max*price_data_low

    group_columns <- c("sim", summarise_by)

    dyn_value <- bind_mse_outputs(model_runs, c("land_caa"), extra_columns) %>%
        as_tibble() %>%
        group_by(across(all_of(c("time", group_columns)))) %>%
        mutate(tot_catch = sum(value)) %>%
        filter(fleet == "Fixed") %>%
        filter_times(time_horizon = time_horizon) %>%
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
        group_by(across(all_of(c("time", group_columns)))) %>%
        summarise(total_value = sum(dyn_price*value)) %>%
        group_by(across(all_of(group_columns))) %>%
        summarise(dyn_annual_value = mean(total_value)) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "dyn_annual_value",
            rel_value = relative,
            grouping = group+columns
        )

    if(!is.null(extra_filter)){
        dyn_value <- dyn_value %>% filter(eval(extra_filter))
    }

    return(
        dyn_value %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(dyn_annual_value, .width=interval_widths, .simple_names=FALSE)
    )

}

#' Compute all Performance Metrics
#' 
#' Wrapper functions to compute all performance metrics (identified below).
#' - average annual catch
#' - total catch
#' - proportion of years catch exceeds a threshold
#' - average ssb
#' - average annual catch variation
#' - average proportion of catch that is "large"
#' - average proportion of population that is "old"
#' - annual relative value
#' - annual dynamic value
#' 
#' The summarised output only includes total catch, average ssb, average 
#' annual catch variation, proportion large catch, proportion old population,
#' and annual dynamic value.
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' 
#' @export performance_metric_summary
#'
#' @example
#'
performance_metric_summary <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    interval_widths,
    time_horizon = c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr")
){
    # Average Catch Across Projection Period
    avg_catch <- average_catch(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Total Catch Across Projection Period
    tot_catch <- total_catch(model_runs, extra_columns, interval_widths, extra_filter=extra_filter, time_horizon=time_horizon, relative=relative, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Prop Years High Catch
    prop_years_high_catch <- prop_years_catch(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, catch_threshold = 30, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Average SSB Across Projection Period
    avg_ssb <- average_ssb(model_runs, extra_columns, dem_params, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Average Annual Catch Variation Across Projection Period
    avg_variation <- average_annual_catch_variation(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Average proportion of catch that is "large"
    avg_catch_lg <- average_proportion_catch_large(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Average proportion of population that is "old"
    avg_pop_old <- average_proportion_biomass_old(model_runs, extra_columns, dem_params, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    # Average annual value
    annual_value <- average_annual_value(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    dynamic_value <- average_annual_dynamic_value(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by) %>% reformat_ggdist_long(n=length(summarise_by))

    perf_data <- bind_rows(avg_catch, avg_ssb, avg_variation, avg_catch_lg, avg_pop_old, dynamic_value) %>%
        mutate(name=factor(
                        name, 
                        levels=c("annual_catch", "total_catch", "num_years", "spbio", "aav", "catch", "bio", "annual_value", "dyn_annual_value"), 
                        labels=c("Annual Catch", "Total Catch", "Years High Catch", "SSB", "Catch AAV", "Large Catch", "Old SSB", "Annual Value", "Dynamic Annual Value")
                    )
        )

    return(afscOM::listN(avg_catch, tot_catch, prop_years_high_catch, avg_ssb, avg_variation, avg_catch_lg, avg_pop_old, annual_value, dynamic_value, perf_data))

}
