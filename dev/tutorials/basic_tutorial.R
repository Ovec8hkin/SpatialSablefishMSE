rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SPoCK)
library(doParallel)

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

# Define a named list of OM and MP objects
# Lists MUST be named, so recommend using afscOM::listN()
om_list <- afscOM::listN(om_rand_recruit, om_bhcyclic_recruit)
hcr_list <- afscOM::listN(mp_f40, mp_20cap)

# Setup MSE options with helper function
mse_options <- setup_mse_options()
# Going to only project 15 years for simplicity
mse_options$n_proj_years <- 15

# Define a number of simulations and simulation seeds
# Parallel processing function can run up to N-1 simulations concurrently
# where N is the number of cores available on your machine.
nsims <- 10
seeds <- sample(1:1e4, nsims, replace=FALSE)

model_runs <- run_mse_multiple(
    om_list = om_list,
    hcr_list = hcr_list,
    mse_options_list = afscOM::listN(mse_options),
    seed_list = seeds,
    diagnostics = FALSE,
    save = FALSE
)

model_runs <- readRDS(file.path(here::here(), "dev", "tutorials", "model_runs_basic.rds"))

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
hcr_colors <- set_hcr_colors2(hcrs)

### Lets Process Output Data -------------------------
ssb_data <- get_ssb_biomass(
    model_runs = model_runs,
    extra_columns = extra_columns,
    dem_params = om_rand_recruit$dem_params,
    hcr_filter = hcrs,
    om_filter = oms
)


# Aggregated SSB data across simulations according to the set interval widths
# and plot as regional depletion and Alaska-wide spawning biomass.
depletion_plots <- plot_depletion(
        ssb_data, 
        v1="hcr", 
        v2="om", 
        v3="region", 
        common_trajectory=common_trajectory, 
        show_est = FALSE, 
        scales="fixed"
    )+
    scale_y_continuous(limits=c(0, 3.5))+
    scale_color_manual(values=c("red", "blue"))+
    labs(title="Regional Depletion (SSB/SSB0)", y="Depletion")

ssb_agg_plots <- plot_ssb(ssb_data, v1="hcr", v2="om", common_trajectory=common_trajectory, show_est=TRUE, scales="free_y", base_hcr = "F45")+
    # scale_y_continuous(limits=c(0, 350))+
    coord_cartesian(ylim=c(0, 350))+
    scale_color_manual(values=c("red", "blue"))+
    facet_wrap(~om, ncol=1)+
    labs(title="Alaska-Wide Spawning Biomass")
    # guides(color="none", shape="none")
    # theme(
    #     strip.background = element_blank(),
    #     strip.text = element_blank()
    # )
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
reg_catch_plots <- plot_landed_catch(catch_data, v1="hcr", v2="om", v3="region", by_fleet=FALSE, common_trajectory = common_trajectory)+
    scale_y_continuous(limits=c(0, 20))+
    scale_color_manual(values=c("red", "blue", "green", "purple"))+
    labs(title="Regional Landed Catch")

catch_agg_plot <- plot_landed_catch(catch_data, v1="hcr", v2="om", by_fleet=TRUE, common_trajectory = common_trajectory)+
    scale_color_manual(values=c("red", "blue", "green", "purple"))+
    facet_wrap(~om, ncol=1)+
    guides(color="none", shape="none")+
    labs(title="Alaska-Wide Landed Catch")+
    # theme(
    #     strip.background = element_blank(),
    #     strip.text = element_blank()
    # )

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
        facet_grid(om~region)+
        custom_theme+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
