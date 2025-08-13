#' Get Spawning Biomass and Total Biomass
#' 
#' Process MSE simulations for spawning biomass,
#' and total stock biomass.
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

    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            filter_hcr_om(hcr_filter, om_filter) %>%
            drop_na() %>%
            # join WAA and maturity-at-age for computing SSB
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
            mutate(
                om = factor(om, levels=om_filter),
                hcr = factor(hcr, levels=hcr_filter)
            )
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
    group_columns <- c("time", "fleet", "region", "sim", "L1", names(extra_columns))
    
    return(
        bind_mse_outputs(model_runs, c("faa", "faa_est"), extra_columns) %>% 
            as_tibble() %>%
            filter(hcr %in% hcr_filter, om %in% om_filter) %>%
            drop_na() %>%
            group_by(across(all_of(group_columns))) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                F = max(value)
            ) %>%
            # ungroup() %>%
            # group_by(across(all_of(group_columns[-2]))) %>%
            # # total F is the sum of fleet-based Fs
            # mutate(
            #     total_F = sum(F)
            # ) %>%
            # ungroup() %>%
            mutate(
                om = factor(om, levels=om_filter),
                hcr = factor(hcr, levels=hcr_filter)
            )
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
    group_columns <- c("time", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            filter(hcr %in% hcr_filter, om %in% om_filter) %>%
            drop_na() %>%
            filter(age == 2) %>%
            group_by(across(all_of(c(group_columns, "region")))) %>%
            summarise(rec=sum(value)) %>%
            mutate(
                om = factor(om, levels=om_filter),
                hcr = factor(hcr, levels=hcr_filter)
            )
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
    group_columns <- c("time", "fleet", "region", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("land_caa"), extra_columns) %>%
            as_tibble() %>%
            filter(hcr %in% hcr_filter, om %in% om_filter) %>%
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
            ungroup() %>%
            mutate(
                om = factor(om, levels=om_filter),
                hcr = factor(hcr, levels=hcr_filter)
            )
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
    cols <- c("time", "sim", "region", "fleet", "value", "L1", names(extra_columns))

    mgmt <- bind_mse_outputs(model_runs, c("abc", "tac", "exp_land"), extra_columns) %>%
                as_tibble() %>%
                filter_hcr_om(hcrs=hcr_filter, oms=om_filter) %>%
                drop_na() %>%
                select(cols) %>%
                pivot_wider(names_from=L1, values_from=value) %>%
                pivot_longer(abc:exp_land, names_to="L1", values_to="value") %>%
                mutate(
                    om = factor(om, levels=om_filter),
                    hcr = factor(hcr, levels=hcr_filter)
                )

    return(mgmt)
}

get_numbers_at_age <- function(model_runs, extra_columns, hcr_filter, om_filter){
    group_columns <- c("time", "class", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa"), extra_columns) %>%
            as_tibble() %>%
            filter(hcr %in% hcr_filter, om %in% om_filter) %>%
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
get_reference_points <- function(model_runs, extra_columns, hcr_filter, om_filter, seed_list){

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

    avg_recruitment <- get_recruits(model_runs, extra_columns, hcr_filter, om_filter) %>%
        group_by(sim, om) %>%
        summarise(rec=mean(rec))

    prop_fs_df <- get_fishing_mortalities(model_runs, extra_columns, hcr_filter, om_filter) %>%
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
        reframe(rps = get_rps(om, hcr, rec, time, c(Fixed, Trawl))) %>%
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

    b40_tseries <- get_rp_timeseries(model_runs, extra_columns, hcr_filter, om_filter, refpt="Bref", spr_target=0.40)

    b40s <- b40_tseries %>% 
        group_by(time, om) %>%
        median_qi(B40, .width=interval_widths)

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

    avg_recruitment <- get_recruits(model_runs, extra_columns, hcr_filter, om_filter) %>% 
        filter(L1 == "naa_est") %>% 
        group_by(sim, om, hcr) %>%
        mutate(avg_rec = unlist(lapply(slider::slide(rec, ~.x, .before=Inf), \(x) mean(x)))) %>%
        arrange(hcr, om, sim)

    prop_fs_df <- get_fishing_mortalities(model_runs, extra_columns, hcr_filter, om_filter) %>%
        filter(L1 != "faa_est") %>%
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
        left_join(avg_recruitment %>% select(-c(L1)), by=c("time", "sim", "om", "hcr")) %>%
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
#' @param n_proj_years number of projection years to compute convergence diagnostics over
#' @param nsims number of simulations over which to calculate convergence diagnostics
#' 
#' @return dataframe of whether EM run for projection year y, and simulation s, converged
#'         according to the above diagnostic checks.
#' 
#' @export check_em_convergence_diagnostics
#' 
check_em_convergence_diagnostics <- function(model_runs, extra_columns, hcr_filter, om_filter, n_proj_years, nsims){

    post_optim_sanity_checks <- function(sd_rep, rep, gradient_tol = 1e-3, se_tol = 100, corr_tol = 0.99) {
        passed_post_sanity_checks <- TRUE
        
        # check likelihoods are all finite and not NA
        if(!all(is.finite(rep$jnLL))) {
            message("Found Inf in joint log-likelihood, model is not converged!")
            passed_post_sanity_checks <- F
        }
        # check maximum absolute gradients
        max_abs_grad_ndx <- which.max(abs(sd_rep$gradient.fixed))
        max_abs_grad <- abs(sd_rep$gradient.fixed)[max_abs_grad_ndx]
        if(gradient_tol < max_abs_grad) {
            message("Parameter: ", names(sd_rep$par.fixed)[max_abs_grad_ndx], " had absolute gradient = ", max_abs_grad,
                    " which was greater than tolerance ", gradient_tol,". This indicates potential non-convergence according to the tolerance.\n")
            passed_post_sanity_checks <- F
        }
        # check hessian
        if(!sd_rep$pdHess) {
            message("Hessian is not positive definite, model is not converged!")
            passed_post_sanity_checks <- F
        }
        # check if standard errors are finite
        if(!all(is.finite(sqrt(diag(sd_rep$cov.fixed))))) {
            message("Found non finite elements in standard errors of parameters, model is not converged!")
            passed_post_sanity_checks <- F
        }
        # check if standard errors are big
        if(max(sqrt(diag(sd_rep$cov.fixed))) > se_tol) {
            message("Parameter: ", names(diag(sd_rep$cov.fixed))[which.max(sqrt(diag(sd_rep$cov.fixed)))], " has a standard error = ",
                    max(sqrt(diag(sd_rep$cov.fixed))), " which was greated than tolerance ", se_tol, ". This indicates potential non-convergence according to the tolerance. \n")
            passed_post_sanity_checks <- F
        }
        # check if correlations are big
        corr_mat <- cov2cor(sd_rep$cov.fixed)
        diag(corr_mat) <- "Same" # set diagonal to "Same" to remove from max calculations
        # reshape to dataframe
        corr_df <- reshape2::melt(corr_mat) %>%
            dplyr::filter(value != 'Same') %>%
            dplyr::mutate(value = as.numeric(value))
        if(max(abs(corr_df$value)) > corr_tol) {
            message("Parameter pairs: ", corr_df$Var1[which.max(abs(corr_df$value))], " and ", corr_df$Var2[which.max(abs(corr_df$value))], " have a correlation of ", max(abs(corr_df$value)), ". This indicates potential non-convergence according to the tolerance.")
            passed_post_sanity_checks <- F
        }
        cat("\n\n");
        if(passed_post_sanity_checks) {
            message("Successfully passed post-optim-sanity checks\n")
        }

        return(passed_post_sanity_checks)

    }

    object_num <- extra_columns %>% rownames_to_column() %>%
        filter_hcr_om(hcr_filter, om_filter) %>% pull(rowname)

    mse_obj <- model_runs[[as.numeric(object_num)]]
    

    # nsims <- length(em_model_obj)/length(em_model_obj[[1]])

    output <- expand.grid(year=1:(n_proj_years+1), sim=1:nsims)
    output$converged <- NA
    i <- 1
    for(s in 1:nsims){
        for(y in 0:n_proj_years){
            em_model_obj <- mse_obj$model_outs[[(y+1)+(s-1)*(n_proj_years+1)]]
            rep <- em_model_obj$rep
            sdrep <- sdreport(em_model_obj)
            converged <- tryCatch(
                post_optim_sanity_checks(sdrep, rep, corr_tol=0.999),
                error = function(e) FALSE
            )
            output[i,]$converged <- converged
            i = i+1
        }
    }
    
    return(output)

}