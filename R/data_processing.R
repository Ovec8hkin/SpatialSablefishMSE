#' Get Spawning Biomass and Total Biomass
#' 
#' Process MSE simulations for spawning biomass,
#' total stock biomass, and depletion. Depletion calculated as the ratio
#' of annual SSB to SSB in timestep 1.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param dem_params demographic parameter matrix
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_ssb_biomass
#'
get_ssb_biomass <- function(model_runs, extra_columns, dem_params, hcr_filter, om_filter){
    process <- function(data){
        data %>%
            left_join(
                melt(dem_params$waa, value.name="weight"),
                by=c("time", "age", "sex", "region")
            ) %>%
            left_join(
                melt(dem_params$mat, value.name="maturity"),
                by=c("time", "age", "sex", "region")
            ) %>%
            drop_na() %>%
            # compute derived quantities
            mutate(
                biomass = value*weight,
                spbio = value*weight*maturity
            ) %>%
            # SSB is females only
            filter(sex == "F") %>%
            # summarise SSB across year and sim 
            group_by(across(all_of(c(group_columns, "region")))) %>%
            summarise(
                spbio=sum(spbio),
                biomass=sum(biomass)
            ) %>%
            ungroup()
    }

    group_columns <- c("time", "sim", "L1", names(extra_columns))

    dem_params$waa <- array(
        dem_params$waa[,,,1], 
        dim=c(dim(dem_params$waa)[1:3], dim(dem_params$waa)[4]+1),
        dimnames=c(dimnames(dem_params$waa)[1:3], list("region"=c("BS", "AI", "WGOA", "CGOA", "EGOA", "Alaska")))
    )

    dem_params$mat <- array(
        dem_params$mat[,,,1], 
        dim=c(dim(dem_params$mat)[1:3], dim(dem_params$mat)[4]+1),
        dimnames=c(dimnames(dem_params$mat)[1:3], list("region"=c("BS", "AI", "WGOA", "CGOA", "EGOA", "Alaska")))
    )

    out <- process_big_outputs(model_runs, c("naa", "naa_est"), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)

    ssb0 <- out %>% filter(time == 1) %>% select(sim, L1, om, hcr, region, ssb0=spbio)

    return(
        out %>% left_join(ssb0, by=c(group_columns[-1], "region")) %>% mutate(dep=spbio/ssb0)
    )
}

#' Get Annual Fishing Mortality
#' 
#' Process MSE simulations for fishing mortality by fleet.
#' Total fishing mortality across fleets is also computed as the
#' sum of fleet-specific fishing mortality.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_fishing_mortalities
#'
get_fishing_mortalities <- function(model_runs, extra_columns, hcr_filter, om_filter){
    process <- function(data){
        data %>%
            group_by(across(all_of(group_columns))) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                F = max(value)
            )
    }
    group_columns <- c("time", "fleet", "region", "sim", "L1", names(extra_columns))
    return(
        process_big_outputs(model_runs, c("faa", "faa_est"), extra_columns, hcr_filter, om_filter, process) %>% 
            format(hcr_filter, om_filter)
    )
}

#' Get Annual Recruits
#' 
#' Process MSE simulations for annual recruits.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_recruits
#'
get_recruits <- function(model_runs, extra_columns, hcr_filter, om_filter){
    process <- function(data){
        data %>%
            filter(age == 2) %>%
            group_by(across(all_of(c(group_columns, "region")))) %>%
            summarise(rec=sum(value))
    }
    group_columns <- c("time", "sim", "L1", names(extra_columns))
    return(
        process_big_outputs(model_runs, c("naa", "naa_est"), extra_columns, hcr_filter, om_filter, process) %>% 
            format(hcr_filter, om_filter)
    )
}

#' Get Landed Catches
#' 
#' Process MSE simulations for landed catches by fleet.
#' Total landed catch across fleets is also computed.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_landed_catch
#'
get_landed_catch <- function(model_runs, extra_columns, hcr_filter, om_filter){
    process <- function(data){
        data %>% group_by(across(all_of(group_columns))) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                catch = sum(value)
            ) %>%
            ungroup() %>%
            group_by(across(all_of(group_columns[-2]))) %>%
            # total F is the sum of fleet-based Fs
            mutate(
                total_catch = sum(catch)
            ) %>%
            ungroup()
    }
    group_columns <- c("time", "fleet", "region", "sim", "L1", names(extra_columns)) 
    return(
        process_big_outputs(model_runs, c("land_caa"), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)
    )
}

#' Get Catch and Population Numbers
#' 
#' Process MSE simulations for landed catches by fleet in numbers of fish.
#' Total landed catch in numbers across fleets is also computed. Numbers
#' caught are calculated as C = N*(1-exp(-F*v)), where N is the number of fish 
#' in the population, F is the fishing mortality rate.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_catch_pop_numbers
#'
get_catch_pop_numbers <- function(model_runs, extra_columns, dem_params, hcr_filter, om_filter){
    process <- function(data){
        data %>% as_tibble() %>%
            pivot_wider(names_from=L1, values_from=value) %>%
            # left_join(
            #     melt(dem_params$sel, value.name="selectivity"),
            #     by=c("time", "age", "sex", "region", "fleet")
            # ) %>%
            mutate(
                true_faa = faa#*selectivity
            ) %>%
            group_by(across(all_of(c(group_columns, "age", "sex")))) %>%
            mutate(faa = sum(true_faa, na.rm=TRUE)) %>%
            filter(is.na(fleet)) %>%
            select(-fleet, -true_faa) %>%
            mutate(
                dead = naa*(1-exp(-faa))
            ) %>%
            group_by(across(all_of(group_columns))) %>%
            summarise(
                total_dead = sum(dead, na.rm=TRUE),
                total_num = sum(naa, na.rm=TRUE)
            )
            # select(-faa, -dead, -naa)
    }
    group_columns <- c("time", "region", "sim", names(extra_columns)) 

    return(
        process_big_outputs(model_runs, c("naa", "faa"), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)
    )
}

#' Get ABC, TAC, and Expected Landings
#' 
#' Process MSE simulations for ABC, TAC, and expected 
#' landings quantities. Historical ABCs, TACs, and landings (e.g., from the
#' conditioning period) are not available.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#' @param spinup_years year at which to begin calculating quantities
#'
#' @export get_management_quantities
#'
get_management_quantities <- function(model_runs, extra_columns, hcr_filter, om_filter, spinup_years=54){
    process <- function(data){
        data %>%
            select(group_columns) %>%
            # filter_times(time_horizon = time_horizon) %>%
            pivot_wider(names_from=L1, values_from=value) %>%
            pivot_longer(abc:exp_land, names_to="L1", values_to="value")
    }
    group_columns <- c("time", "sim", "region", "fleet", "value", "L1", names(extra_columns))
    return(
        process_big_outputs(model_runs, c("abc", "tac", "exp_land"), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)
    )
}

#' Get Timeseries of Dynamic Economic Value
#' 
#' Process MSE simulations for dynamic economic value (Goethel et al. 2025).
#' Economic value calculated as a linear adjustment in price by size grade
#' between a max price achieved at low total landings and a min price achieved
#' at high total landings. Value is calculated only for the Fixed gear fleet. 
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_dynamic_economic_value
#'
get_dynamic_economic_value <- function(model_runs, extra_columns, hcr_filter, om_filter){
    process <- function(data){
        data %>%
            group_by(across(all_of(c("time", group_columns)))) %>%
            mutate(tot_catch = sum(value)) %>%
            # filter(fleet == "Fixed") %>%
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
            group_by(across(all_of(c("time", "region", "fleet", group_columns)))) %>%
            summarise(total_value = sum(dyn_price*value))
    }

    group_columns <- c("sim", "om", "hcr")

    price_age_f_low <- c(0.597895623, 1.320303448, 1.320303448, 1.856562267, 2.610111345, 2.610111345, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875)
    price_age_m_low <- c(0.597895623, 0.597895623, 1.320303448, 1.320303448, 1.856562267, 1.856562267, 1.856562267, 1.856562267, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345)
    price_data_low <- matrix(c(price_age_f_low, price_age_m_low), nrow=length(price_age_f_low), ncol=2)
    dimnames(price_data_low) <- list("age"=2:31, "sex"=c("F", "M"))

    price_age_f_max <- c(7.917460094, 8.40756497, 8.40756497, 9.944657109, 11.46480347, 11.46480347, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658)
    price_age_m_max <- c(7.917460094, 7.917460094, 8.40756497, 8.40756497, 9.944657109, 9.944657109, 9.944657109, 9.944657109, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347)
    price_data_max <- matrix(c(price_age_f_max, price_age_m_max), nrow=length(price_age_f_max), ncol=2)
    dimnames(price_data_max) <- list("age"=2:31, "sex"=c("F", "M"))
    
    # max_price <- max(c(price_age_f_max, price_age_m_max))
    # price_data_low <- price_data_low/price_data_max
    # price_data_max <- price_data_max/max_price
    # price_data_low <- price_data_max*price_data_low

    return(
        process_big_outputs(model_runs, c("land_caa"), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)
    )

}

#' Get Numbers-at-Price Grade Timeseries
#' 
#' Process MSE simulations for number of fish in each price grade
#' both in catch and in the population. Price grades are defined by
#' age groups as follows:
#' Grade 1/2: <3 yo
#' Grade 2/3: <5 yo
#' Grade 3/4: <7 yo
#' Grade 4/5: <9 yo
#' Grade 5/7: <15 yo
#' Grade 7+: 15+ yo
#'
#' Note that for CAA, units are in weight of catch in each price grade,
#' while for NAA, units are number of fish in each price grade.
#' 
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_numbers_at_age
#'
get_numbers_at_age <- function(model_runs, extra_columns, hcr_filter, om_filter, pop_or_catch="pop"){
    process <- function(data){
        data %>%
            mutate(
                class = factor(
                    case_when(age < 3 ~ "1/2", age < 5 ~ "2/3", age < 7 ~ "3/4", age < 9 ~ "4/5", age < 15 ~ "5/7", age > 14 ~ "7+"), 
                    levels=c("1/2", "2/3", "3/4", "4/5", "5/7", "7+"), 
                    labels=c("Grade 1/2 (1-2yo)", "Grade 2/3 (3-4yo)", "Grade 3/4 (5-6yo)", "Grade 4/5 (7-8yo)", "Grade 5/7 (9-14yo)", "Grade 7+ (15+yo)")
                ),
                L1 = factor(L1, levels=c("caa", "naa"))
            ) %>%
            group_by(across(all_of(group_columns))) %>%
            summarise(value=sum(value))
    }

    v = ifelse(pop_or_catch == "pop", "naa", "caa")

    group_columns <- c("time", "class", "region", "sim", "L1", names(extra_columns))
    return(
        process_big_outputs(model_runs, c(v), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)
            
    )
}

#' Get Average Population Age Timeseries
#' 
#' Process MSE simulations for average age of an individual fish
#' in the population. This is distinct from the performance metric
#' of the same name in that include both male and female fish in the
#' calculation of average age.
#' 
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
#' @export get_average_age
#'
get_average_age <- function(model_runs, extra_columns, hcr_filter, om_filter){
    process <- function(data){
        data %>% filter_times(time_horizon=time_horizon) %>%
            ungroup() %>%
            group_by(time, sim, age, hcr, om, region) %>%
            summarise(value = sum(value)) %>%
            group_by(time, sim, hcr, om, region) %>%
            summarise(
                avg_age = compute_average_age(value, 2:31)
            )
    }

    return(
        process_big_outputs(model_runs, c("naa"), extra_columns, hcr_filter, om_filter, process) %>%
            format(hcr_filter, om_filter)
    )
}

#' Get Catch or Numbers-at-Age by Age Groups
#' 
#' Process MSE simulation results for -at-age by
#' specified age groups. 
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param q either "caa" (for catch-at-age) or "naa" (for numbers-at-age)
#' @param age_groups ages that define age groups (3 groups required)
#' @param group_names names for each age group
#' @param group_abbs abbreviated names for each age group
#' @param summarise whether to summarise across simualtions or not
#' @param make_segments whether to generate data.frame of segment for use in plotting 
#'
#'
get_atage_groups <- function(model_runs, extra_columns, hcr_filter, om_filter, q, age_groups, group_names, group_abbs, spinup_years=64, summarise=FALSE, make_segments=FALSE){
    data <- bind_mse_outputs(model_runs, c(q), extra_columns) %>%
            as_tibble() %>%
            filter(hcr %in% hcr_filter, om %in% om_filter) %>%
            mutate(
                class = factor(
                    case_when(age < age_groups[1] ~ group_names[1], age < age_groups[2] ~ group_names[2], TRUE ~ group_names[3]), 
                    levels=group_names, 
                    labels=group_abbs
                ),
                L1 = factor(L1, levels=c("caa", "naa"), labels=c("Catch-at-Age", "Numbers-at-Age"))
            ) %>%
            group_by(time, class, sim, L1, hcr, om) %>%
            summarise(value=sum(value)) %>%
            filter(time > spinup_years) %>%
            select(time, class, sim, hcr, om, value) %>%
            pivot_wider(names_from="hcr", values_from="value") %>%
            group_by(time, sim, om) %>%
            mutate(across(3:(ncol(.)-3), ~ ./sum(.))) %>% 
            pivot_longer(6:(ncol(.)), names_to="hcr", values_to="value") %>%
            ungroup() %>%
            pivot_wider(names_from="class", values_from="value")

    if(summarise){
        data <- data %>% group_by(time, hcr, om) %>%
            filter(time > spinup_years) %>%
            summarise(across((ncol(.)-2-3):(ncol(.)-3), ~ mean(.)))
    }

    if(make_segments){
        segments <- data %>% as_tibble() %>% 
            select(2:ncol(.)) %>% 
            rename(x=3, y=4, z=5) %>%
            group_by(om, hcr) %>%
            mutate(
                xend = lead(x, 1),
                yend = lead(y, 1),
                zend = lead(z, 1)
            ) %>%
            ungroup() %>%
            arrange(om, hcr) %>%
            drop_na()
        
        return(afscOM::listN(data, segments))
    }else{
        return(data)
    }
}

#' Get Reference Points from MSE Simulations
#' 
#' Derive fishing mortality and biomass reference points
#' from completed MSE simulations.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be appended to the resultant data frame
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#' @param seed_list simulation seeds used in `model_runs`
#'
get_reference_points <- function(model_runs, extra_columns, hcr_filter, om_filter, seed_list, year=55){

    om_names <- om_filter
    hcr_names <- hcr_filter

    get_rps <- function(om_name, hcr_name, recruitment, year, prop_fs){
        om <- om_list[which(om_names == om_name)]
        hcr <- hcr_list[which(hcr_names == hcr_name)]

        om <- om[[1]]
        hcr <- hcr[[1]]

        joint_selret <- calculate_joint_selret(
            sel=om$dem_params$sel[year,,,,,drop=FALSE],
            ret=om$dem_params$ret[year,,,,,drop=FALSE],
            prop_fs = prop_fs
        )
        ref_pts <- calculate_ref_points(
            30,
            om$dem_params$mort[year,,1,1],
            om$dem_params$mat[year,,1,1],
            om$dem_params$waa[year,,1,1],
            joint_selret$sel[,,1,,drop=FALSE],
            joint_selret$ret[,,1,,drop=FALSE],
            recruitment/2,
            spr_target = hcr$ref_points$spr_target
        )
        return(c(ref_pts$Fref, ref_pts$Fmax, ref_pts$Bref, ref_pts$B0))
    }

    r <- get_recruits(model_runs, extra_columns, hcr_filter, om_filter)
    groups <- colnames(r)
    groups <- groups[!groups %in% c("region", "rec")]

    avg_recruitment <- r %>%
        filter(L1 == "naa") %>%
        group_by(across(all_of(groups))) %>%
        summarise(rec=sum(rec)) %>%
        group_by(sim, om) %>%
        summarise(rec=mean(rec))

    fs <- get_fishing_mortalities(model_runs, extra_columns, hcr_filter, om_filter)

    prop_fs_df <- fs %>%
        filter(L1 == "faa_est") %>%
        pivot_wider(names_from="fleet", values_from="F") %>%
        mutate(total_F = Fixed + Trawl) %>%
        pivot_longer(Fixed:Trawl, names_to="fleet", values_to="F") %>%
        group_by(time, sim, om, hcr, fleet) %>%
        mutate(
            prop_f = F/total_F
        ) %>%
        select(time, sim, fleet, om, hcr, prop_f) %>%
        distinct() %>%
        pivot_wider(names_from = "fleet", values_from="prop_f") %>%
        group_by(sim, om, hcr) %>%
        summarise(Fixed = mean(Fixed), Trawl=mean(Trawl))


    ref_pts_df <- prop_fs_df %>% 
        left_join(avg_recruitment, by=c("sim", "om")) %>%
        group_by(sim, om, hcr) %>%
        reframe(rps = get_rps(om, hcr, rec, year, c(Fixed, Trawl))) %>%
        mutate(rp_name = rep(c("Fref", "Fmax", "Bref", "B0"), length(hcr_filter)*length(om_filter)*length(seed_list))) %>%
        pivot_wider(names_from="rp_name", values_from="rps") %>%
        group_by(om, hcr) %>%
        median_qi(Fref, Fmax, Bref, B0, .width=interval_widths, .simple_names=TRUE) %>%
        reformat_ggdist_long(n=2) %>% 
        filter(.width == 0.5) %>%
        select(om, hcr, name, median) %>%
        pivot_wider(names_from=name, values_from=median) %>%
        arrange(hcr)


    return(ref_pts_df)

}

#' Get Timeseries of B40 from MSE Simulations
#' 
#' Derive timeseries of B40 from completed MSE
#' simulations.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#' 
#' @return timeseries of B40 reference point summarised across simulations
#'
get_b40_timeseries <- function(model_runs, extra_columns, hcr_filter, om_filter){

    b40_tseries <- get_rp_timeseries(model_runs, extra_columns, hcr_filter, om_filter, ref_pt="Bref", spr_target=0.40)

    b40s <- b40_tseries %>% 
        group_by(time, om) %>%
        median_qi(rp, .width=interval_widths)

    return(b40s)

}


#' Get Timeseries of Refernce Point from MSE Simulations
#' 
#' Derive timeseries of a given reference point from completed MSE
#' simulations.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#' @param ref_pt the reference point to return (one of "Fmax", "Fref", "Bref", or "B0")
#' @param spr_target target SPR reference level to compute timeseries of (if "hcr", will calculate
#' uniquely for each HCR based on the the SPR target reference point defined in the HCR object)
#' 
#' @return timeseries of reference point value
#'
get_rp_timeseries <- function(model_runs, extra_columns, hcr_filter, om_filter, ref_pt, spr_target="hcr"){

    get_rps <- function(om_name, hcr_name, recruitment, year, prop_fs){

        om <- om_list[which(om_filter == om_name)]
        hcr <- hcr_list[which(hcr_filter == hcr_name)]

        om <- om[[1]]
        hcr <- hcr[[1]]

        joint_selret <- calculate_joint_selret(
            sel=om$dem_params$sel[year,,,,,drop=FALSE],
            ret=om$dem_params$ret[year,,,,,drop=FALSE],
            prop_fs = prop_fs
        )

        ref_pts <- calculate_ref_points(
            30,
            om$dem_params$mort[year,,1,1],
            om$dem_params$mat[year,,1,1],
            om$dem_params$waa[year,,1,1],
            joint_selret$sel[,,1,,drop=FALSE],
            joint_selret$ret[,,1,,drop=FALSE],
            recruitment/2,
            spr_target = ifelse(spr_target=="hcr", hcr$ref_points$spr_target, spr_target)
        )
        return(ref_pts[[ref_pt]])
    }

    r <- get_recruits(model_runs, extra_columns, hcr_filter, om_filter)
    groups <- colnames(r)
    groups <- groups[!groups %in% c("region", "rec")]

    avg_recruitment <- r %>% 
        filter(L1 == "naa") %>% 
        group_by(across(all_of(groups))) %>%
        summarise(rec=sum(rec)) %>%
        group_by(sim, om, hcr) %>%
        mutate(avg_rec = unlist(lapply(slider::slide(rec, ~.x, .before=Inf), \(x) mean(x)))) %>%
        arrange(hcr, om, sim)

    fs <- get_fishing_mortalities(model_runs, extra_columns, hcr_filter, om_filter)

    prop_fs_df <- fs %>%
        filter(L1 == "faa_est") %>%
        pivot_wider(names_from="fleet", values_from="F") %>%
        mutate(total_F = Fixed + Trawl) %>%
        pivot_longer(Fixed:Trawl, names_to="fleet", values_to="F") %>%
        group_by(time, sim, om, hcr, fleet) %>%
        mutate(
            prop_f = F/total_F
        ) %>%
        select(time, sim, fleet, om, hcr, prop_f) %>%
        distinct() %>%
        pivot_wider(names_from = "fleet", values_from="prop_f") %>%
        group_by(time, sim, om, hcr) %>%
        summarise(Fixed = mean(Fixed), Trawl=mean(Trawl))

    rps <- prop_fs_df %>% 
        left_join(avg_recruitment, by=c("time", "sim", "om", "hcr")) %>%
        rowwise() %>%
        mutate(rp = get_rps(om, hcr, avg_rec, time, c(Fixed, Trawl))) #%>% 
        # group_by(time, om) %>%
        # median_qi(B40, .width=interval_widths)

    return(rps)

}

#' Get F and SSB data for Phaseplane Plot
#' 
#' Process MSE simulation data for total F
#' and SSB by year across simulations, OMs,
#' and management procedures.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param dem_params list of demographic parameters values
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
get_phaseplane_data <- function(model_runs, extra_columns, dem_params, hcr_filter, om_filter){
    return(
        get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
            ungroup() %>%
            select(-L1) %>%
            left_join(
                get_fishing_mortalities(model_runs, extra_columns, hcr_filter, om_filter) %>% filter(L1 == "faa", fleet == "Fixed") %>% select(time, sim, om, hcr, total_F),
                by = c("time", "sim", "hcr", "om"),
            )
    )
}

#' Get HCR recommended F and SSB data for Phaseplane Plot
#' 
#' Process MSE simulation data for recommended 
#' F (from HCR) and SSB by year across simulations, 
#' OMs, and management procedures.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param dem_params list of demographic parameters values
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
get_hcrphase_data <- function(model_runs, extra_columns, dem_params, hcr_filter, om_filter){
    return(
        get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
            # SSB is females only
            # filter(sex == "F", L1 == "naa") %>%
            # summarise SSB across year and sim 
            # group_by(time, hcr, sim, L1) %>%
            # summarise(spbio=sum(spbio)) %>%
            ungroup() %>%
            select(-L1) %>%
            left_join(
                bind_mse_outputs(model_runs, "hcr_f", extra_columns) %>% as_tibble() %>% filter(hcr %in% hcr_filter, om %in% om_filter),
                by = c("time", "sim", "hcr", "om"),
            ) %>%
            select(-c(Var2, Var3, region, L1))
    )
}

#' Get Catch and SSB data for Phaseplane Plot
#' 
#' Process MSE simulation data for landed catch
#' and SSB by year across simulations, OMs, and 
#' management procedures.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param dem_params list of demographic parameters values
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#'
get_phaseplane_catch_data <- function(model_runs, extra_columns, dem_params, hcr_filter, om_filter){
    return(
        get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>%
            ungroup() %>%
            select(-L1) %>%
            left_join(
                get_landed_catch(model_runs, extra_columns, hcr_filter, om_filter) %>% filter(L1 == "land_caa", fleet == "Fixed") %>% select(time, sim, om, hcr, total_catch),
                by = c("time", "sim", "hcr", "om"),
            )
    )
}

#' Check EM Convergence Diagnostics
#' 
#' Checks a returned EM object for convergence based on the following properties:
#' - No infinite log-likelihood components
#' - Maximum gradient < 1e-3
#' - PD hessian
#' - No inifinte SEs of parameters
#' - No SEs > 100
#' - No parameter correlations > 0.999
#' 
#' @param model_runs list of completed MSE model runs with `diagnostics=TRUE`
#' @param extra_columns additional columns to associated model with HCR and OM
#' @param hcr_filter vector of HCR names to process (must match names in `extra_columns`)
#' @param om_filter vector of OM names to process (must match names in `extra_columns`)
#' 
#' @return dataframe of whether EM run for projection year y, and simulation s, converged
#'         according to the above diagnostic checks.
#' 
#' @export get_convergence_data
#' 
get_convergence_data <- function(model_runs, extra_columns, hcr_filter, om_filter){

    fs <- list.files(file.path(here::here(), "data", "active"), full.names = TRUE)
    if(!is.null(hcr_filter))
        fs <- unlist(sapply(hcr_filter, \(x) fs[grepl(paste0(sub("|", "\\|", sub("/", "", sub(" ", "_", tolower(x))), fixed=TRUE), "_\\d+"), fs)]))

    o <- tibble::tibble()

    k <- 1
    for(i in seq_along(fs)){
        x <- fs[i]
        print(x)
        m <- readRDS(x)
        mse <- m$mse_objects
        model_run <- list(mse[[length(mse)]])[[1]]
        extra_columns <- expand.grid(om=model_run$om$name, hcr=model_run$mp$name)

        if(!(model_run$om$name %in% om_filter)){
            next;
        }

        seeds <- m$seeds

        convergence <- model_run$converged

        fullconvergence <- reshape2::melt(convergence, value.name="converged") %>% 
            as_tibble() %>%
            rename(sim_year="Var1", sim="Var2") %>%
            mutate(
                om=extra_columns$om, 
                hcr=extra_columns$hcr,
                sim = factor(sim, labels=seeds)
            )
        
        o <- bind_rows(o, fullconvergence)
        print(o %>% nrow())

    }
    return(o)
}

#' Calculation Bias Between OM and EM Outputs
#' 
#' Calculates relative between quantities from the OM and estimates from the 
#' EM. Because the EM is spatially aggregated, biases are calculated at the
#' Alaska-wide level, rather than at the regional level. Bias is calculated
#' as (estimate - true)/true, where the true value is the sum of the regional
#' values from the OM.
#' 
#' This function has only been tested on SSB and recruitment, 04/21/2026.
#' 
#' @param data tibble of output data with both OM and EM outputs
#' @param om_em_cols values of L1 column in data indicating OM and EM outputs (e.g., "naa" and "naa_est")
#' @param value_col name of column in data containing values to compare between OM and EM
#' @param time_horizon vector of time steps to include in bias calculation
#' 
#' @return tibble of bias between OM and EM outputs
#' 
#' @export get_estimation_bias
get_estimation_bias <- function(data, om_em_cols, value_col, time_horizon){
    data %>% select(time, sim, L1, om, Recruitment, Movement, hcr, region, value_col) %>%
        filter_times(time_horizon) %>%
        mutate(region = ifelse(is.na(region), "Alaska", as.character(region))) %>%
        pivot_wider(names_from=L1, values_from=value_col) %>%
        group_by(time, sim, om, Recruitment, Movement, hcr) %>%
        mutate(true = sum(eval(rlang::parse_expr(om_em_cols[1])), na.rm=TRUE)) %>%
        rename(est = om_em_cols[2]) %>%
        filter(region == "Alaska") %>%
        ungroup() %>%
        mutate(
            bias = (est - true)/true
        ) %>%
        select(-om_em_cols[1])
}
