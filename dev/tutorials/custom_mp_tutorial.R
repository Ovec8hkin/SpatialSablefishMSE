library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SPoCK)
library(doParallel)
library(patchwork)

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

devtools::load_all(afscOM_dir)

# Load the SpatialSablefishMSE functions
lapply(list.files("R", full.names = TRUE), source)

# Load predefined OM and HCRs
source("dev/oms.R")   
source("dev/hcrs.R")

om_list <- listN(om_rand_recruit)

### Lets define a new harvest control rule ----------

# Start from a base MP object so we dont miss any required fields
mp_new <- setup_mp_options()

# Define reference points
target_spr <- c(0.50, 0.40) # c(Ftar, Btar) --> c(F50, B40)
mp_new$ref_points$spr_target <- target_spr

calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,1],
    mat = dp_y$mat[,,1,1],
    waa = dp_y$waa[,,1,1],
    sel =  joint_selret$sel[,,1,1,drop=FALSE],
    ret = joint_selret$ret[,,1,1,drop=FALSE],
    avg_rec = mean(apply(hist_recruits/2, 1, sum)),
    spr_target = target_spr
)

# Define the harvest control rule

#' All HCR functions MUST take the following arguments:
#' - ref_pts: a list of reference points (e.g. Fmax, Bref, Fref, B0)
#' - naa: the current numbers at age array ([1, nages, nsexes, nregions])
#' - dem_params: demographic parameters (e.g. waa, mat; [1, nages, nsexes, nregions, nfleets])
#' - avgrec: average recruitment (not used in this example, but required)
#
#' Can probably get around the required arguments by using `...` in the function definition,
#' but this is not recommended as it can lead to confusion about what arguments are required.
new_hcr <- function(ref_pts, naa, dem_params, avgrec, cutoff_age=1, lrp=0.05){
    nages <- afscOM::get_model_dimensions(naa)$nages
    ssb <- apply(naa[,cutoff_age:nages,1,,drop=FALSE]*dem_params$waa[,cutoff_age:nages,1,,drop=FALSE]*dem_params$mat[,cutoff_age:nages,1,,drop=FALSE], 1, sum)
    return(
        threshold_f(
            x = ssb / ref_pts$Bref, 
            f_min = 0, 
            f_max = ref_pts$Fref, 
            lrp = lrp, 
            urp = 1
        )
    )
}

mp_new$hcr
mp_new$hcr$func <- new_hcr

# Set the value of the additional HCR function parameter
# `cutoff_age`.
mp_new$hcr$extra_pars <- list(
    cutoff_age = 3, # only calculate SSB for ages 4+,
    lrp = 0.20   # LRP 0.20 of Bref
)

# Set additional HCR options like stability constraints and 
# harvest caps.
mp_new$hcr$extra_options <- list(
    # Define stability constraints c(%down, %up)
    max_stability = array(c(1.0, 0.05), dim=c(1, 2)),
    # Define harvest cap (1000s mt)
    harvest_cap = 20
)
# HCR function output is in fishing mortality units
mp_new$hcr$units <- "F"

# Define aspects of management policy
mp_new$management
mp_new$management$abc_tac_regflt_reduction <- array(
    # Give BS trawl an extra 10% TAC at the expense of AI trawl
    matrix(
        c(1, 1.10, 1, 0.90, 1, 1, 1, 1, 1, 1),
        nrow=5,ncol=2, byrow=TRUE
    ),
    dim = c(5, 2),
    dimnames = list(
        "region" = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
        "fleet" = c("Fixed", "Trawl")
    )
)
# Assume full TAC utilization
mp_new$management$regflt_tac_utilization <- full_utilization

# Don't forget a regional apportionment function

# All apportionment functions MUST take the following arguments:
# - survey_obs: a survey observation object (e.g. from a survey model)
# - y: the year for which to apportion
#
# Can probably get around the required arguments by using `...` in the function definition,
# but this is not recommended as it can lead to confusion about what arguments are required.
#
# Apportionment functions should return a vector/array of dimensions [1, nregions]
apportionment <- function(survey_obs, y, window_size=5){
    rpw <- survey_obs$rpws[1:y,,,,1]
    roll <- apply(rpw, 2, \(x) zoo::rollmean(x, k=window_size))
    app <- t(apply(roll, 1, \(x) x/sum(x)))
    return(app[nrow(app),])
}

mp_new$apportionment$func <- apportionment
mp_new$apportionment$pars <- list(
    window_size = 3 # Rolling mean window size of 3 years
)

mp_new$name <- "Hybrid HCR"

hcr_list <- listN(mp_f40, mp_new)

### Setup the everything else we need for the MSE -----------
mse_options <- setup_mse_options()
mse_options$n_proj_years <- 15

nsims <- 10
seeds <- sample(1:1e4, nsims, replace=FALSE)

model_runs <- run_mse_multiple(
    om_list = om_list,
    hcr_list = hcr_list,
    mse_options_list = listN(mse_options),
    seed_list = seeds,
    diagnostics = FALSE,
    save = FALSE
)  

# Going to read them from the file instead of running the MSE for time purposes
model_runs <- readRDS(file.path(here::here(), "dev", "tutorials", "model_runs_custom_mp.rds"))

oms <- unlist(lapply(om_list, \(x) x$name)) # Names of the OMs you want to use in data processing
hcrs <- unlist(lapply(hcr_list, \(x) x$name)) # Names of the HCRs you want to use in data processing

# Grid of OM/HCR combinations for model_runs
extra_columns <- expand.grid(om=oms, hcr=hcrs)

# Quantiles for distributions
interval_widths <- c(0.50, 0.80) 
# How many years of trajectory are shared across OMs/HCRs
common_trajectory <- mse_options$n_spinup_years 
# Time horizon over which to calculate metrics
time_horizon <- c(mse_options$n_spinup_years+1, mse_options$n_spinup_years+mse_options$n_proj_years) 

# HCR colors for plots
hcr_colors <- c("F40"="red", "Hyrbid HCR"="blue")

### Lets Process Output Data -------------------------

# Get SSB data for all OMs, HCRs, and regions across all simulations
ssb_data <- get_ssb_biomass(
    model_runs = model_runs,
    extra_columns = extra_columns,
    dem_params = om_rand_recruit$dem_params,
    hcr_filter = hcrs,
    om_filter = oms
)

# Aggregated SSB data across simulations according to the set interval widths
# and plot as regional depletion and Alaska-wide spawning biomass.
depletion_plots <- plot_depletion(ssb_data, v1="hcr", v2="om", v3="region", common_trajectory=common_trajectory, show_est = FALSE, scales="fixed")+
    scale_y_continuous(limits=c(0, 3.5))+
    scale_color_manual(values=c("red", "blue"))+
    labs(title="Regional Depletion (SSB/SSB0)", y="Depletion")

ssb_agg_plots <- plot_ssb(ssb_data, v1="hcr", v2="om", common_trajectory=common_trajectory, scales="free_y")+
    scale_y_continuous(limits=c(0, 350))+
    scale_color_manual(values=c("red", "blue"))+
    facet_wrap(~om, ncol=1)+
    labs(title="Alaska-Wide Spawning Biomass")+
    guides(color="none", shape="none")+
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )
# Combine the plots into a single figure
ssb_agg_plots + depletion_plots + 
    plot_layout(nrow=1, guides="collect", widths=c(0.5, 1)) & 
    theme(legend.position = "bottom")



# Get landings data for all OMs, HCRs, and regions across all simulations
catch_data <- get_landed_catch(
    model_runs = model_runs,
    extra_columns = extra_columns,
    hcr_filter = hcrs,
    om_filter = oms
)

# Aggregate landings data across simulations according to the set interval widths
# and plot as regional landings and Alaska-wide landings.
reg_catch_plots <- plot_landed_catch(catch_data, v1="hcr", v2="om", v3="region", by_fleet=TRUE, common_trajectory = 20)+
    scale_y_continuous(limits=c(0, 20))+
    scale_color_manual(values=c("red", "blue", "green", "purple"))+
    labs(title="Regional Landed Catch")
    # ggh4x::facet_nested(om +fleet~ region)
catch_agg_plot <- plot_landed_catch(catch_data, v1="hcr", v2="om", by_fleet=TRUE, common_trajectory = common_trajectory)+
    scale_color_manual(values=c("red", "blue", "green", "purple"))+
    facet_wrap(~om, ncol=1)+
    guides(color="none", shape="none")+
    labs(title="Alaska-Wide Landed Catch")+
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )

catch_agg_plot + reg_catch_plots + 
    plot_layout(nrow=1, guides="collect", axes = "collect", widths=c(0.5, 1)) & 
    theme(legend.position = "bottom")



# Get ABC, TAC, and expected landings data for all OMs and HCRs
abc_tac_land <- get_management_quantities(
    model_runs = model_runs,
    extra_columns = extra_columns,
    hcr_filter = hcrs,
    om_filter = oms,
    spinup_years = mse_options$n_spinup_years
)

abc_tac_land_reg <- abc_tac_land %>%
    mutate(
        L1 = factor(L1, levels=c("abc", "tac", "exp_land"), labels=c("ABC", "TAC", "Landed Catch")),
    ) %>%
    group_by(time, region, om, hcr, fleet, L1) %>%
    median_qi(value, .width=interval_widths) %>%
    filter(.width == 0.50)

regional_abctac_plots <- ggplot(abc_tac_land_reg, aes(x=time, y=value, color=hcr, linetype=L1))+
    geom_line()+
    ggh4x::facet_nested(om + fleet ~ region, scales="free_y")+
    scale_color_manual(values=c("red", "blue"))+
    custom_theme+
    labs(title="Regional ABC, TAC, and Landed Catch")+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 15)),
            scale_y_continuous(limits=c(0, 5)),
            scale_y_continuous(limits=c(0, 15)),
            scale_y_continuous(limits=c(0, 8))
        )
    )

abc_tac_land_agg <- abc_tac_land %>%
    mutate(
        L1 = factor(L1, levels=c("abc", "tac", "exp_land"), labels=c("ABC", "TAC", "Landed Catch")),
    ) %>%
    group_by(time, om, hcr, fleet, L1, sim) %>%
    summarise(value=sum(value)) %>%
    group_by(time, om, hcr, fleet, L1) %>%
    median_qi(value, .width=interval_widths) %>%
    filter(.width == 0.50)

agg_abctac_plot <- ggplot(abc_tac_land_agg, aes(x=time, y=value, color=hcr, linetype=L1))+
    geom_line()+
    facet_wrap(om ~ fleet, scales="free_y", ncol=1)+
    labs(title="Alaska-Wide ABC, TAC, and Landed Catch")+
    scale_color_manual(values=c("red", "blue"))+
    custom_theme+
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 35)),
            scale_y_continuous(limits=c(0, 10)),
            scale_y_continuous(limits=c(0, 40)),
            scale_y_continuous(limits=c(0, 15))
        )
    )

agg_abctac_plot + regional_abctac_plots + 
    plot_layout(nrow=1, guides="collect", axes = "collect", widths=c(0.5, 1)) & 
    theme(legend.position = "bottom")


### Calculate Performance Metrics -----------------
performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    # Compute metrics across OMs, HCRs, and regions
    summarise_by=c("om", "hcr", "region"),
    hcr_filter=hcrs,
    om_filter=oms,
    # Compute average catch, average SSB, average female age, and average dynamic value metrics
    # More metrics available in the roxygen documentation
    metric_list = c("avg_catch", "avg_ssb", "avg_age", "dynamic_value")
)

# Human readable metric names
publication_metrics = c("Annual Catch", "SSB", "Average Age", "Dynamic Annual Value")
publication_metrics2 = c("Annual\nCatch", "SSB", "Average\nAge", "Dynamic\nValue")

perf_data <- performance_metrics$perf_data %>% 
    filter_hcr_om(hcrs, oms) %>%
    mutate(
        hcr = factor(hcr, levels=hcrs),
        om = factor(om, labels=oms),
        name = factor(name, levels=publication_metrics, labels=publication_metrics2)
    ) %>%
    filter(.width == 0.5) %>%
    group_by(om, region, name) %>%
    # scale relative to the maximum value for each metric to obtain 0-1 range
    scale_and_rank("median")

perf_data %>% select(om, hcr, region, name, scaled) %>%
    ggplot(aes(x=name, y=scaled, color=hcr, fill=hcr, group=hcr))+
        geom_point(size=3)+
        geom_line()+
        geom_polygon(alpha=0)+
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))+
        ggiraphExtra::coord_radar()+
        facet_wrap(~region)+
        custom_theme+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
