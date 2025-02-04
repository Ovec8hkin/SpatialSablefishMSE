#' Compute Average Catch across projection period
#' 
#' Compute average catch (median and CIs) per year across all years and simulation seeds, 
#' for each combination of operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_catch
#'
#' @example
#'
average_catch <- function(
    model_runs, 
    extra_columns, 
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA),  
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){

    group_columns <- c("sim", summarise_by)
    
    avg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_hcr_om(hcr_filter, om_filter) %>%
        filter_times(time_horizon) %>%
        group_by(across(all_of(c("time", group_columns)))) %>%
        summarise(annual_catch = sum(value)) %>%
        round_to_zero("annual_catch") %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "annual_catch",
            rel_value = relative,
            grouping=group_columns
        )
    
    if(!is.null(extra_filter)){
        avg_catch <- avg_catch %>% filter(eval(extra_filter))
    }

    out <- avg_catch
    if(summary_out){
        out <- avg_catch %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(annual_catch, .width=interval_widths, .simple_names=FALSE)
    }
            
    return(
        out 
    )
}

#' Compute Total Catch over projection period
#' 
#' Compute total catch (median and CIs) over all years and simulation seeds, 
#' for each combination of operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export total_catch
#'
#' @example
#'
total_catch <- function(
    model_runs, 
    extra_columns, 
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    
    group_columns <- c("sim", summarise_by)
    tot_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_hcr_om(hcr_filter, om_filter) %>%
        filter_times(time_horizon = time_horizon) %>%
        group_by(across(all_of(group_columns))) %>%
        summarise(total_catch = sum(value)) %>%
        round_to_zero("total_catch") %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "total_catch",
            rel_value = relative,
            grouping=group_columns
        )

    if(!is.null(extra_filter)){
        tot_catch <- tot_catch %>% filter(eval(extra_filter))
    }

    out <- tot_catch
    if(summary_out){
        out <- tot_catch %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(total_catch, .width=interval_widths, .simple_names=FALSE)
    }
            
    return(out)
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
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export prop_years_catch
#'
#' @example
#'
prop_years_catch <- function(
    model_runs, 
    extra_columns, 
    hcr_filter,
    om_filter,
    catch_threshold, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE
){
    
    group_columns <- c("sim", summarise_by)
    catch_years <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_hcr_om(hcr_filter, om_filter) %>%
        filter_times(time_horizon = time_horizon) %>%
        group_by(across(all_of(c("time", group_columns)))) %>%
        summarise(
            total_catch = sum(value),
        ) %>%
        round_to_zero("total_catch") %>%
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

    out <- catch_years
    if(summary_out){
        out <- catch_years %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(num_years, .width=interval_widths, .simple_names=FALSE)
    }
            
    return(out)
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
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_ssb
#'
#' @example
#'
average_ssb <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE
){
    
    group_columns <- c("sim", summarise_by)

    avg_ssb <- get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
            ungroup() %>%
            filter(L1 != "naa_est") %>%
            select(-c("L1", "biomass")) %>%
            filter_times(time_horizon=time_horizon) %>%
            round_to_zero("spbio") %>%
            relativize_performance(
                rel_column = "hcr",
                value_column = "spbio",
                rel_value = relative,
                grouping = group_columns
            )
            
    if(!is.null(extra_filter)){
        avg_ssb <- avg_ssb %>% filter(eval(extra_filter))
    }

    out <- avg_ssb
    if(summary_out){
        out <- avg_ssb %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(spbio, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
}

#' Compute proportion of years SSB is below a threshold level across projection period
#' 
#' Compute average proportion of years where SSB is below 0.4*SSB0
#' (median and CIs) across all years and simulation seeds, 
#' for each combination of operating models and management procedures. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param time_horizon periof of years over which to calculate metric
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export prop_low_biomass
#'
#' @example
#'
prop_low_biomass <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(55, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    all_ssb <- get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
                        ungroup() %>%
                        filter(L1 != "naa_est") %>%
                        select(-c("L1", "biomass")) %>%
                        round_to_zero("spbio")
    
    group_columns <- c("sim", summarise_by)

    #threshold <- all_ssb %>% filter(time == 1) %>% pull(spbio) %>% min
    threshold <- 0.35*299.901

    prop_years_low_biomass <- all_ssb %>% 
        filter_times(time_horizon=time_horizon) %>%
        mutate(
            low_bio = spbio <= threshold
        ) %>%
        group_by(sim, om, hcr) %>%
        summarise(
            total_low_bio = sum(low_bio),
            prop_years = total_low_bio/(time_horizon[2]-time_horizon[1])
        ) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "prop_years",
            rel_value = relative,
            grouping = group_columns
        )
    
    if(!is.null(extra_filter)){
        prop_years_low_biomass <- prop_years_low_biomass %>% filter(eval(extra_filter))
    }

    out <- prop_years_low_biomass
    if(summary_out){
        out <- prop_years_low_biomass %>% 
            group_by(across(all_of(summarise_by))) %>%
            median_qi(prop_years, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
}

#' Compute number of years requied for SSB to crash
#' 
#' Compute average number of years required for SSB to decline below
#' 0.2*SSB0 (median and CIs) for each combination of operating 
#' models and management procedures.
#' 
#' Scenarios (combinations of OM-MP-seed) where the crash threshold is not
#' met, are removed from the dataset before computing median and CIs. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export biomass_crash_time
#'
#' @example
#'
biomass_crash_time <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(55, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    all_ssb <- get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
                        ungroup() %>%
                        filter(L1 != "naa_est") %>%
                        select(-L1) %>%
                        round_to_zero("spbio")

    #threshold <- all_ssb %>% filter(time == 1) %>% pull(spbio) %>% min
    threshold <- 0.50*0.35*299.901

    group_columns <- c("sim", summarise_by)

    crash_time <- all_ssb %>% filter_times(time_horizon=c(time_horizon[1], time_horizon[1]+20)) %>%
        filter(spbio <= threshold) %>%
        group_by(sim, om, hcr) %>%
        summarise(crash_time=min(time)-time_horizon[1]) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "crash_time",
            rel_value = NULL,
            grouping = group_columns
        )
    
    if(!is.null(extra_filter)){
        crash_time <- crash_time %>% filter(eval(extra_filter))
    }

    out <- crash_time
    if(summary_out){
        out <- crash_time %>% 
            group_by(across(all_of(summarise_by))) %>%
            median_qi(crash_time, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
}

#' Compute number of years requied for SSB to recover
#' 
#' Compute average number of years required for SSB to recover above
#' 0.4*SSB0 (median and CIs) following a known SSB crash
#' for each combination of operating models and management procedures. 
#' 
#' "Recovery" begins exactly 20 years after the start of the simulation
#' to correspond with the "Immediate Crash Recruitment" scenario. 
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export biomass_recovery_time
#'
#' @example
#'
biomass_recovery_time <- function(
    model_runs, 
    extra_columns, 
    dem_params, 
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(55, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    all_ssb <- get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
                        ungroup() %>%
                        filter(L1 != "naa_est") %>%
                        select(-L1) %>%
                        round_to_zero("spbio")

    # threshold <- all_ssb %>% filter(time == 1) %>% pull(spbio) %>% min
    threshold <- 0.35*299.901

    group_columns <- c("sim", summarise_by)

    recovery_time <- all_ssb %>% filter_times(time_horizon=c(time_horizon[1]+20, time_horizon[2])) %>%
        filter(spbio >= threshold) %>%
        group_by(sim, om, hcr) %>%
        summarise(recovery_time=min(time)-(time_horizon[1]+20)) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "recovery_time",
            rel_value = NULL,
            grouping = group_columns
        )
    
    if(!is.null(extra_filter)){
        recovery_time <- recovery_time %>% filter(eval(extra_filter))
    }

    out <- recovery_time
    if(summary_out){
        out <- recovery_time %>% 
            group_by(across(all_of(summarise_by))) %>%
            median_qi(recovery_time, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
}

#' Compute Average Age across projection period
#' 
#' Compute average age (median and CIs) per year across all years and 
#' simulation seeds, for each combination of operating models and management 
#' procedures. Average age is calculated by age of individuals and is not
#' weighted by biomass or maturity.
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param dem_params demographic parameters matrices from OM
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_age
#'
#' @example
#'
average_age <- function(
    model_runs, 
    extra_columns,
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    
    group_columns <- c("sim", summarise_by)

    avg_age <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
            as_tibble() %>%
            ungroup() %>%
            filter_hcr_om(hcr_filter, om_filter) %>%
            filter_times(time_horizon=time_horizon) %>%
            group_by(time, age, sim, om, hcr) %>%
            mutate(value = sum(value)) %>%
            # round_to_zero("value") %>%
            filter(sex == "F") %>%
            ungroup() %>%
            group_by(time, sim, hcr, om) %>%
            summarise(
                avg_age = compute_average_age(value, 2:31)
            ) %>%
            relativize_performance(
                rel_column = "hcr",
                value_column = "avg_age",
                rel_value = relative,
                grouping = group_columns
            )
            
    if(!is.null(extra_filter)){
        avg_age <- avg_age %>% filter(eval(extra_filter))
    }

    out <- avg_age
    if(summary_out){
        out <- avg_age %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(avg_age, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
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
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_abi
#'
#' @example
#'
average_abi <- function(
    model_runs, 
    extra_columns,
    ref_naa,
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    
    group_columns <- c("sim", summarise_by)

    avg_abi <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
            as_tibble() %>%
            ungroup() %>%
            filter_hcr_om(hcr_filter, om_filter) %>%
            filter_times(time_horizon=time_horizon) %>%
            group_by(time, age, sim, om, hcr) %>%
            mutate(value = sum(value)) %>%
            # round_to_zero("value") %>%
            filter(sex == "F") %>%
            ungroup() %>%
            group_by(time, sim, hcr, om) %>%
            summarise(
                avg_abi = abi(value, ref_naa)
            ) %>%
            relativize_performance(
                rel_column = "hcr",
                value_column = "avg_abi",
                rel_value = relative,
                grouping = group_columns
            )
            
    if(!is.null(extra_filter)){
        avg_abi <- avg_abi %>% filter(eval(extra_filter))
    }

    out <- avg_abi
    if(summary_out){
        avg_abi <- avg_abi %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(avg_abi, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
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
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_annual_catch_variation
#'
#' @example
#'
average_annual_catch_variation <- function(
    model_runs, 
    extra_columns,
    hcr_filter,
    om_filter,
    interval_widths=c(0.50, 0.80),
    time_horizon = c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    
    group_columns <- c("sim", summarise_by)
    avg_catch_var <- get_landed_catch(model_runs, extra_columns, hcr_filter, om_filter) %>%
        filter_times(time_horizon = time_horizon) %>%
        round_to_zero("total_catch") %>%
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

    out <- avg_catch_var
    if(summary_out){
        out <- avg_catch_var %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(aav, .width=interval_widths, .simple_names=FALSE)
    }

    return(out)
}

#' Compute average proportion of catch that is "Large" across projection period
#' 
#' Compute average proportion of catch that is "large" (median and CIs) 
#' per year across all years and simulation seeds, for each combination 
#' of operating models and management procedures. "Large" fish are those 
#' >9yo (the approximate age corresponding to the highest processor grade
#' for sablefish in Alaska).
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_proportion_catch_large
#'
#' @example
#'
average_proportion_catch_large <- function(
    model_runs, 
    extra_columns,
    hcr_filter,
    om_filter, 
    interval_widths=c(0.50, 0.80),
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    
    group_columns <- c("sim", summarise_by)
    prop_lg_catch <- bind_mse_outputs(model_runs, "caa", extra_columns) %>%
        as_tibble() %>%
        filter_hcr_om(hcr_filter, om_filter) %>%
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
        round_to_zero("total_catch") %>%
        mutate(across(Large:Small, ~ ./total_catch)) %>%
        round_to_zero("Large") %>%
        round_to_zero("Medium") %>%
        round_to_zero("Small") %>%
        select(-total_catch) %>%
        ungroup() %>%
        pivot_longer(Large:Small, names_to="size_group", values_to="catch") %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "catch",
            rel_value = relative,
            grouping = c("size_group", group_columns)
        ) %>%
        filter(size_group == "Large")

    if(!is.null(extra_filter)){
        prop_lg_catch <- prop_lg_catch %>% filter(eval(extra_filter))
    }

    out <- prop_lg_catch
    if(summary_out){
        out <- prop_lg_catch %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(catch, .width=interval_widths, .simple_names=FALSE)
    }
    
    return(out)
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
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' 
#' @export average_proportion_biomass_old
#'
#' @example
#'
average_proportion_biomass_old <- function(
    model_runs, 
    extra_columns, 
    dem_params,
    hcr_filter,
    om_filter, 
    interval_widths=c(0.50, 0.80), 
    time_horizon=c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE
){
    
    group_columns <- c("sim", summarise_by)
    prop_old_biomass <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
        as_tibble() %>%
        filter_hcr_om(hcr_filter, om_filter) %>%
        filter_times(time_horizon = time_horizon) %>%
        # join WAA and maturity-at-age for computing SSB
        left_join(
            melt(dem_params$waa, value.name="weight"), 
            by=c("time", "age", "sex")
        ) %>%
        mutate(bio = value*weight) %>%
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
        round_to_zero("total_bio") %>%
        mutate(across(Adult:Young, ~ ./total_bio)) %>%
        round_to_zero("Young") %>%
        round_to_zero("Adult") %>%
        round_to_zero("Old") %>%
        select(-total_bio) %>%
        ungroup() %>%
        pivot_longer(Adult:Young, names_to="age_group", values_to="bio") %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "bio",
            rel_value = relative,
            grouping = c("age_group", group_columns)
        ) %>%
        filter(age_group == "Old")

    if(!is.null(extra_filter)){
        prop_old_biomass <- prop_old_biomass %>% filter(eval(extra_filter))
    }

    out <- prop_old_biomass
    if(summary_out){
        out <- prop_old_biomass %>%
            group_by(across(all_of(summarise_by))) %>%
            median_qi(bio, .width=interval_widths, .simple_names=FALSE)
    }
    
    return(out)
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
#' @param summary_out whether to output data summarised by `ggdist` or full data
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
            grouping = group_columns
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
#' Wrapper function to compute all performance metrics (identified below).
#' - average annual catch
#' - total catch
#' - proportion of years catch exceeds a threshold
#' - average ssb
#' - proportion of years SSB is below a threshold
#' - average age of the population
#' - average ABI of the population (Griffiths et al. 2023)
#' - average annual catch variation
#' - average proportion of catch that is "large"
#' - average proportion of population that is "old"
#' - annual relative value
#' - annual dynamic value
#' - number of years required for population to delcine below 0.2B0
#' - number of years required for population to exceed below 0.4B0
#'
#' @param model_runs list of completed MSE simulations runs
#' @param extra_columns data.frame specifying names for OM and HCR to attach
#' to each model_run (see `bind_mse_outputs` for more details)
#' @param dem_params demographic parameters matrices from OM
#' @param ref_naa reference age structure for ABI computation
#' @param hcr_filter vector of HCR names to calculate metric over
#' @param om_filter vector of OM names to calculate metric over
#' @param interval_widths confidence intevrals to compute
#' @param extra_filter an additional set of filters to apply before computing 
#' medians and confidence intervals
#' @param relative a management procedure to compute metric relative to
#' @param summarise_by vector of columns to summarise metric by
#' @param summary_out whether to output data summarised by `ggdist` or full data
#' @param metric_list vector of names of performance metrics to compute
#'      - `avg_catch` -> average annual catch
#'      - `tot_catch` -> total catch
#'      - `prop_years_high_catch` -> proportion of years catch exceeds a threshold
#'      - `avg_ssb` -> average ssb
#'      - `prop_years_lowssb` -> proportion of years SSB is below a threshold
#'      - `avg_age` -> average age of the population
#'      - `avg_abi` -> average ABI of the population (Griffiths et al. 2023)
#'      - `avg_variation` -> average annual catch variation
#'      - `avg_catch_lg` -> average proportion of catch that is "large"
#'      - `avg_pop_old` -> average proportion of population that is "old"
#'      - `annual_value` -> annual relative value
#'      - `dynamic_value` -> annual dynamic value
#'      - `crash_time` -> number of years required for population to delcine below 0.2B0
#'      - `recovery_time` -> number of years required for population to exceed below 0.4B0
#'      - `all` -> compute all of the above metrics
#' 
#' @return list of dataframes. One element for each individual metric, and one element with
#' all metrics
#' 
#' @export performance_metric_summary
#'
#' @example
#'
performance_metric_summary <- function(
    model_runs, 
    extra_columns, 
    dem_params,
    ref_naa, 
    hcr_filter,
    om_filter,
    interval_widths,
    time_horizon = c(65, NA), 
    extra_filter=NULL, 
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out=TRUE,
    metric_list="all"
){

    n <- ifelse(summary_out, length(summarise_by), 0)

    # Average Catch Across Projection Period
    if(any(c("avg_catch", "all") %in% metric_list))
        avg_catch <- average_catch(model_runs, extra_columns, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Total Catch Across Projection Period
    if(any(c("tot_catch", "all") %in% metric_list))
        tot_catch <- total_catch(model_runs, extra_columns, hcr_filter, om_filter, interval_widths, extra_filter=extra_filter, time_horizon=time_horizon, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Prop Years High Catch
    if(any(c("prop_years_high_catch", "all") %in% metric_list))
        prop_years_high_catch <- prop_years_catch(model_runs, extra_columns, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, catch_threshold = 30, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average SSB Across Projection Period
    if(any(c("avg_ssb", "all") %in% metric_list))
        avg_ssb <- average_ssb(model_runs, extra_columns, dem_params, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average SSB Across Projection Period
    if(any(c("prop_years_lowssb", "all") %in% metric_list))
        prop_years_lowssb <- prop_low_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average SSB Across Projection Period
    if(any(c("avg_age", "all") %in% metric_list))
        avg_age <- average_age(model_runs, extra_columns, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average SSB Across Projection Period
    if(any(c("avg_abi", "all") %in% metric_list))
        avg_abi <- average_abi(model_runs, extra_columns, ref_naa, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average Annual Catch Variation Across Projection Period
    if(any(c("avg_variation", "all") %in% metric_list))
        avg_variation <- average_annual_catch_variation(model_runs, extra_columns, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average proportion of catch that is "large"
    if(any(c("avg_catch_lg", "all") %in% metric_list))
        avg_catch_lg <- average_proportion_catch_large(model_runs, extra_columns, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average proportion of population that is "old"
    if(any(c("avg_pop_old", "all") %in% metric_list))
        avg_pop_old <- average_proportion_biomass_old(model_runs, extra_columns, dem_params, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Average annual value
    if(any(c("annual_value", "all") %in% metric_list))
        annual_value <- average_annual_value(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    # Dynamic annual value
    if(any(c("dynamic_value", "all") %in% metric_list))
        dynamic_value <- average_annual_dynamic_value(model_runs, extra_columns, interval_widths, time_horizon=time_horizon, relative=relative, extra_filter=extra_filter, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    if(any(c("crash_time", "all") %in% metric_list))
        crash_time <- biomass_crash_time(model_runs, extra_columns, dem_params, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    if(any(c("recovery_time", "all") %in% metric_list))
        recovery_time <- biomass_recovery_time(model_runs, extra_columns, dem_params, hcr_filter, om_filter, interval_widths, time_horizon=time_horizon, extra_filter=extra_filter, relative=relative, summarise_by = summarise_by, summary_out=summary_out) %>% reformat_ggdist_long(n=n)

    
    metric_key <- c(
        "Annual Catch"="annual_catch", 
        "Total Catch"="total_catch", 
        "Years of High Catch"="num_years", 
        "SSB"="spbio", 
        "Proportion of Years with Low SSB"="prop_years", 
        "Catch AAV"="aav", 
        "Proportion Large Catch"="catch", 
        "Proportion Old Biomass"="bio", 
        "Annual Value"="annual_value", 
        "Dynamic Annual Value"="dyn_annual_value", 
        "Average ABI"="avg_abi", 
        "Average Age"="avg_age",
        "Crash Time"="crash_time",
        "Recovery Time"="recovery_time"
    )  

    perf_data <- NULL  

    if(summary_out){
        perf_data <- bind_rows(mget(metric_list)) %>%
            mutate(name=factor(
                            name, 
                            levels=metric_key, 
                            labels=names(metric_key)
                        )
            )
    }else{
        dfs <- mget(metric_list)
        have_time <- all(unlist(lapply(dfs, \(x) "time" %in% colnames(x))))
        if(all(have_time)){
            perf_data <- plyr::join_all(dfs, by=c("time", "sim", "om", "hcr"))
        }else{
            dfs <- lapply(dfs[1:length(dfs)], function(x){
                x %>% group_by(sim, om, hcr) %>%
                    summarise(across(ncol(.)-3, \(y) median(y, na.rm=TRUE)))
            })
            perf_data <- plyr::join_all(dfs, by=c("sim", "om", "hcr"))
        }
    }

    return(append(mget(metric_list), list(perf_data=perf_data)))

}
