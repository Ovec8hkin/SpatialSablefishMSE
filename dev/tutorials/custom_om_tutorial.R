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

hcr_list <- afscOM::listN(mp_f40)

### Let's define a new OM with no movement

# Set everything to default (which includes a movement matrix)
om_nomove <- om_base

# Get OM dimensions from the base model (assuming we aren't changing
# the dimensionality of the OM at all).
model_dims <- afscOM::get_model_dimensions(om_nomove$dem_params$sel)

# Construct a simple movement matrix with 1s on the diagonal and 0s off the diagonal
nomovement_matrix <- diag(1, nrow=model_dims$nregions, ncol=model_dims$nregions)

# Expand simple movement matrix to the correct dimensions of [nregions, nregions, 
# nyears, nages, nsexes]. We are just going to borrow the dimensions and dimnames 
# from the existing movement matrix since they are the same in this case.
nomovement <- array(
    nomovement_matrix, 
    dim=dim(om_nomove$dem_params$movement), 
    dimnames=dimnames(om_nomove$dem_params$movement)
)

# overwrite the old movement matrix
om_nomove$dem_params$movement[,,55:200,,] <- nomovement[,,55:200,,]

om_nomove$recruitment <- om_rand_recruit$recruitment


om_nomove$name <- "No Movement"





### Let's define a new OM with spatially varying selectivity for the fixed gear fleet
om_spatsel <- sable_om

# Going to make a new array of selectivity-at-age for each region
# for the fixed gear fleet. These selectivity patterns aren't
# biologically realistic, but should serve as an illustrative example
# of how to modify parameters across space.
new_selex <- array(
    c(
        # knife-edged at age-6 in region 1 (BS)
        rep(0, 4), rep(1, 26), 
        # linear increase and then plateau at age 9 (AI)
        0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, rep(1, 23), 
        # semi-logistic (WGOA)
        0, 0.1, 0.35, 0.6, 0.8, 0.9, 0.95, 0.98, 0.99, rep(1, 21), 
        # semi dome-shaped (CGOA)
        0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, rep(1, 15), 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 
        # fully selected (EGOA)
        rep(1, 30) 
    ),
    dim = c(30, 5),
    dimnames = dimnames(om_spatsel$dem_params$sel)[c("age", "region")]
)

# Generate a new parameter matrix of the correct dimensions using the
# afscOM helper function.
selex_dp <- afscOM::generate_param_matrix(
    # array defining our new selectivities
    new_selex,
    # dimension names of the output (should match with all of the
    # other demographic matrices, and must match the dimnames
    # of `new_selex` in the "age" and "region" dimension)
    dimension_names = c(dimnames(om_spatsel$dem_params$sel)[c("time", "age", "sex", "region")], list("fleet"="Fixed")),
    # Dimensions over which `new_selex` varies (row, col)
    by = c("age", "region"),
    # Include the fifth fleet dimension in the output for simplicity
    include_fleet_dim = TRUE
)

# Replace the seectivity of th fixed gear fleet with our new
# parameter matrix
om_spatsel$dem_params$sel[55:200,,,,1] <- selex_dp[55:200,,,,,drop=FALSE]

om_spatsel$recruitment$func <- beverton_holt
om_spatsel$recruitment$pars <- list(
    h = 0.8,
    R0 = 15,
    S0 = 300,
    sigR = 1.20
)
om_spatsel$recruitment$apportionment <- list(
    func = resample_recruit_apportionment,
    pars = list(
        hist_recruits = hist_recruits
    )
)

om_spatsel$name <- "Spatial Selectivity"




om_list <- afscOM::listN(om_rand_recruit, om_nomove, om_spatsel)

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

model_runs <- readRDS(file.path(here::here(), "dev", "tutorials", "model_runs_custom_om.rds"))

# Going to read them from the file instead of running the MSE for time purposes
# model_runs <- readRDS(file.path(here::here(), "dev", "tutorials", "model_runs_custom_mp.rds"))

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
depletion_plots <- plot_depletion(ssb_data, v1="om", v2="hcr", v3="region", common_trajectory=common_trajectory, show_est = FALSE, scales="fixed")+
    scale_y_continuous(limits=c(0, 3.5))+
    scale_color_manual(values=c("red", "blue", "purple"))+
    labs(title="Regional Depletion (SSB/SSB0)", y="Depletion")

ssb_agg_plots <- plot_ssb(ssb_data, v1="om", v2="hcr", common_trajectory=common_trajectory, scales="free_y")+
    scale_y_continuous(limits=c(0, 350))+
    scale_color_manual(values=c("red", "blue", "purple"))+
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
reg_catch_plots <- plot_landed_catch(catch_data, v1="om", v2="hcr", v3="region", by_fleet=TRUE, common_trajectory = common_trajectory)+
    scale_y_continuous(limits=c(0, 40))+
    scale_color_manual(values=c("red", "blue", "purple"))+
    labs(title="Regional Landed Catch")
    # ggh4x::facet_nested(om +fleet~ region)
catch_agg_plot <- plot_landed_catch(catch_data, v1="om", v2="hcr", by_fleet=TRUE, common_trajectory = common_trajectory)+
    scale_color_manual(values=c("red", "blue", "purple"))+
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

regional_abctac_plots <- ggplot(abc_tac_land_reg, aes(x=time, y=value, color=om, linetype=L1))+
    geom_line()+
    ggh4x::facet_nested(hcr + fleet ~ region, scales="free_y")+
    scale_color_manual(values=c("red", "blue", "purple"))+
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

agg_abctac_plot <- ggplot(abc_tac_land_agg, aes(x=time, y=value, color=om, linetype=L1))+
    geom_line()+
    facet_wrap(hcr ~ fleet, scales="free_y", ncol=1)+
    labs(title="Alaska-Wide ABC, TAC, and Landed Catch")+
    scale_color_manual(values=c("red", "blue", "purple"))+
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
    group_by(hcr, region, name) %>%
    # scale relative to the maximum value for each metric to obtain 0-1 range
    scale_and_rank("median")

perf_data %>% select(om, hcr, region, name, scaled) %>%
    ggplot(aes(x=name, y=scaled, color=om, fill=om, group=om))+
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









# Plot Recruitment
recruit_data <- get_recruits(model_runs, extra_columns, hcr_filter=hcrs, om_filter=oms)
reg_rec_plots <- plot_recruitment(recruit_data, v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)+
    scale_color_manual(values=c("red", "blue", "green", "purple"))

recruit_agg_plot <- plot_recruitment(recruit_data, v1="hcr", v2="om", common_trajectory = common_trajectory)+
    scale_color_manual(values=c("red", "blue", "green", "purple"))+
    scale_y_continuous("Recruits (millions)", limits=c(0, 80))+
    facet_wrap(~om, ncol=1)+
    labs(title="Alaska-Wide Recruitment")+
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )

recruit_agg_plot + reg_rec_plots + 
    plot_layout(nrow=1, guides="collect", axes = "collect", widths=c(0.5, 1)) & 
    theme(legend.position = "bottom")
