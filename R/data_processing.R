#' Get Spawning Biomass and Total Biomass
#' 
#' Process MSE simulations for spawning biomass,
#' and total stock biomass.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#' @param dem_params demographic parameter matrix
#'
#' @export get_ssb_biomass
#'
#' @example \dontrun{
#'      dem_params <- om$dem_params
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_ssb_biomass(model_runs, extra_columns, dem_params)
#' }
#'
get_ssb_biomass <- function(model_runs, extra_columns, dem_params){
    group_columns <- c("time", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            # join WAA and maturity-at-age for computing SSB
            left_join(
                melt(dem_params$waa, value.name="weight"), 
                by=c("time", "age", "sex")
            ) %>%
            left_join(
                melt(dem_params$mat, value.name="maturity"), 
                by=c("time", "age", "sex")
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
            group_by(across(all_of(group_columns))) %>%
            summarise(spbio=sum(spbio))
    )
}

#' Get Annual Fishing Mortality
#' 
#' Process MSE simulations for fishing mortality by fleet.
#' Total fishing mortality across fleets is alsoc computed.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_fishing_mortalities
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_fishing_mortalities(model_runs, extra_columns)
#' }
#'
get_fishing_mortalities <- function(model_runs, extra_columns){
    group_columns <- c("time", "fleet", "sim", "L1", names(extra_columns))
    
    return(
        bind_mse_outputs(model_runs, c("faa", "faa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            group_by(across(all_of(group_columns))) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                F = max(value)
            ) %>%
            ungroup() %>%
            group_by(across(all_of(group_columns[-2]))) %>%
            # total F is the sum of fleet-based Fs
            mutate(
                total_F = sum(F)
            ) %>%
            ungroup()
    )
}

#' Get Annual Recruits
#' 
#' Process MSE simulations for annual recruits.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_recruits
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_recruits(model_runs, extra_columns)
#' }
#'
get_recruits <- function(model_runs, extra_columns){
    group_columns <- c("time", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            filter(age == 2) %>%
            group_by(across(all_of(group_columns))) %>%
            summarise(rec=sum(value))
    )
}

#' Get Landed Catches
#' 
#' Process MSE simulations for landed catches by fleet.
#' Total landed catch across fleets is alsoc computed.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_landed_catch
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_landed_catch(model_runs, extra_columns)
#' }
#'
get_landed_catch <- function(model_runs, extra_columns){
    group_columns <- c("time", "fleet", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("land_caa"), extra_columns) %>%
            as_tibble() %>%
            drop_na() %>%
            group_by(across(all_of(group_columns))) %>%
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
    )
}

#' Get ABC, TAC, and Expected Landings
#' 
#' Process MSE simulations for ABC, TAC, and expected 
#' landings quantities. Historical ABCs, TACs, and landings
#' appendded to begining of output data.frame, and are from
#' Goethel et al. 2023.
#'
#' @param model_runs list of completed MSE simulation objects
#' @param extra_columns additional columns to append to output
#' @param spinup_years year at which to begin calculating quantities
#'
#' @export get_management_quantities
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_management_quantities(model_runs, extra_columns)
#' }
#'
get_management_quantities <- function(model_runs, extra_columns, spinup_years=64){
    cols <- c("time", "sim", "value", "L1", names(extra_columns))

    hist_abcs <- c(44200, 37100, 33400, 28800, 25200, 25000, 28800, 25300, 19600, 17200, 16800, 15900, 17200, 16900, 17300, 20900, 23000, 21000, 21000, 20100, 18000, 16100, 15200, 16000, 17200, 16200, 13700, 13700, 11800, 13100, 15000, 15100, 22000, 29600, 34500, 40500)
    hist_tacs <- c(18000, 19300, 17300, 14500, 14800, 13500, 21400, 27700, 36400, 32200, 33200, 28800, 25200, 25000, 28800, 25300, 19400, 16800, 16800, 15400, 17200, 16900, 17300, 20900, 22600, 21000, 20700, 20100, 18000, 16100, 15200, 16000, 17200, 16200, 13700, 13700, 11800, 13100, 15000, 15100, 18300, 26100, 34500, 39600)
    hist_land <- c(10400, 12600, 12000, 11800, 14100, 14500, 28900, 35200, 38400, 34800, 30200, 26400, 23900, 25400, 23600, 20700, 17400, 14600, 13900, 13600, 15600, 14100, 14700, 16400, 17500, 16600, 15600, 16000, 14600, 13100, 11900, 13000, 13900, 13600, 11500, 10900, 10200, 12300, 14200, 16600, 19000, 21300, 26900, 20400)
    historical_management <- data.frame(
        time = (1980:2023)-1960+1,
        abc = c(rep(NA, length(hist_tacs) - length(hist_abcs)), hist_abcs/1000),
        tac = c(rep(NA, length(hist_tacs) - length(hist_tacs)), hist_tacs/1000),
        exp_land = c(rep(NA, length(hist_tacs) - length(hist_land)), hist_land/1000)
    ) %>% as_tibble() %>%
    filter(time <= spinup_years) %>%
    mutate(attainment = exp_land/tac) %>%
    pivot_longer(abc:attainment, names_to="L1", values_to="value")

    mgmt <- bind_mse_outputs(model_runs, c("abc", "tac", "exp_land"), extra_columns) %>%
                as_tibble() %>%
                drop_na() %>%
                select(cols) %>%
                pivot_wider(names_from=L1, values_from=value) %>%
                mutate(attainment = exp_land/tac) %>%
                pivot_longer(abc:attainment, names_to="L1", values_to="value")

    return(
        bind_rows(
            cross_join(mgmt %>% distinct(across(all_of(c("sim", colnames(extra_columns))))), historical_management) %>% filter(time < spinup_years+1),
            mgmt
        )
    )
}

get_numbers_at_age <- function(model_runs, extra_columns){
    group_columns <- c("time", "class", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa"), extra_columns) %>%
            as_tibble() %>%
            mutate(
                class = factor(
                    case_when(age < 3 ~ "1/2", age < 5 ~ "2/3", age < 7 ~ "3/4", age < 9 ~ "4/5", age < 15 ~ "5/7", age > 14 ~ "7+"), 
                    levels=c("1/2", "2/3", "3/4", "4/5", "5/7", "7+"), 
                    labels=c("Grade 1/2 (1-2yo)", "Grade 2/3 (3-4yo)", "Grade 3/4 (5-6yo)", "Grade 4/5 (7-8yo)", "Grade 5/7 (9-14yo)", "Grade 7+ (15+yo)")
                ),
                L1 = factor(L1, levels=c("caa", "naa"), labels=c("Catch-at-Age", "Numbers-at-Age"))
            ) %>%
            group_by(across(all_of(group_columns))) %>%
            summarise(value=sum(value))
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
#' @export get_atage_groups
#'
#' @example
#'
get_atage_groups <- function(model_runs, extra_columns, q, age_groups, group_names, group_abbs, spinup_years=64, summarise=FALSE, make_segments=FALSE){
    data <- bind_mse_outputs(model_runs, c(q), extra_columns) %>%
            as_tibble() %>%
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
#' Note that this function hasn't been tested when multiple
#' MSE simulations are present in the `model_runs` list.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param om_list list of OM objects
#' @param om_names vector of OM names
#'
#' @export get_reference_points
#'
#' @example
#'
get_reference_points <- function(model_runs, extra_columns, om_list, hcr_list){

    om_names <- unlist(lapply(om_list, \(x) x$name))
    hcr_names <- unlist(lapply(hcr_list, \(x) x$name))

    get_rps <- function(om_name, hcr_name, recruitment, prop_fs){
        om <- om_list[[which(om_names == om_name)]]
        hcr <- hcr_list[[which(hcr_names == hcr_name)]]
        year <- 64
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

    avg_recruitment <- get_recruits(model_runs, extra_columns) %>%
        group_by(sim, om) %>%
        summarise(rec=mean(rec))

    prop_fs_df <- get_fishing_mortalities(model_runs, extra_columns) %>%
        filter(L1 != "faa_est") %>%
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
        reframe(rps = get_rps(om, hcr, rec, c(Fixed, Trawl))) %>%
        mutate(rp_name = rep(c("Fref", "Fmax", "Bref", "B0"), length(hcr_list)*length(om_list)*length(seed_list))) %>%
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

get_phaseplane_data <- function(model_runs, extra_columns, dem_params){
    return(
        get_ssb_biomass(model_runs, extra_columns, dem_params) %>%
            ungroup() %>%
            select(-L1) %>%
            left_join(
                get_fishing_mortalities(model_runs, extra_columns) %>% filter(L1 == "faa", fleet == "Fixed") %>% select(time, sim, om, hcr, total_F),
                by = c("time", "sim", "hcr", "om"),
            )
    )
}

get_hcrphase_data <- function(model_runs, extra_columns, dem_params){
    return(
        get_ssb_biomass(model_runs, extra_columns, dem_params) %>%
            # SSB is females only
            # filter(sex == "F", L1 == "naa") %>%
            # summarise SSB across year and sim 
            # group_by(time, hcr, sim, L1) %>%
            # summarise(spbio=sum(spbio)) %>%
            ungroup() %>%
            select(-L1) %>%
            left_join(
                bind_mse_outputs(model_runs, "hcr_f", extra_columns),
                by = c("time", "sim", "hcr", "om"),
            ) %>%
            select(-c(Var2, Var3, region, L1))
    )
}
