rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

devtools::load_all(afscOM_dir)

lapply(list.files("R", full.names = TRUE), source)

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons
nyears <- 110

sable_om <- readRDS("data/sablefish_om_big.RDS") # Read this saved OM from a file
sable_om$model_options$fleet_apportionment <- matrix(c(0.80, 0.20), nrow=nrow(sable_om$model_options$fleet_apportionment), ncol=2, byrow=TRUE)
sable_om$model_options$obs_pars$catch_cv <- c(1e-5, 1e-5, 1e-5, 1e-5)

# Source all available OM and HCR objects
source("dev/oms.R")
source("dev/hcrs.R")

#' 3. Run the closed-loop MSE simulation
#' A single MSE simulation can be run using the `run_mse(...)`
#' function, while multiple MSE simulations can be run (serially)
#' using the `run_mse_multiple(...)` function.
#' 
#' It is recommended to always use `run_mse_multiple(...)` even
#' when only a single MSE simulation is required.
#' 

set.seed(895)
nsims <- 200
seed_list <- sample(1:(100*nsims), nsims)  # Draw 10 random seeds

mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 54
mse_options$recruitment_start_year <- 54
mse_options$n_proj_years <- 75

mse_options_list <- listN(mse_options)


om_list <- listN(om_rand_recruit, om_bh_recruit, om_bhcyclic_recruit, om_immcrash_recruit)
hcr_list <- listN(
    mp_f40, mp_f50,
    mp_5perc, mp_10perc,
    mp_15cap, mp_25cap,
    mp_f50chr,
    mp_pfmc4010, mp_bcsable,
    mp_f00chr
)


tic()
model_runs <- run_mse_multiple(
    om_list, 
    hcr_list, 
    seed_list,
    nyears=100,
    mse_options_list=mse_options_list,
    diagnostics = FALSE,
    save=TRUE
)
toc()


 # Data Processinf
filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures")
width_small <- 12
height_small <- 8

publication_hcrs <- c("F40", "F50", "F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "25k Harvest Cap", "Constant F50", "PFMC 40-10", "British Columbia", "No Fishing")
publication_oms <- c("Random Recruitment", "Beverton-Holt Cyclic Recruitment", "Immediate Crash Recruitment")
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years with Low SSB")

mse_runs <- get_saved_model_runs2(om_order=publication_oms, hcr_order=publication_hcrs)
model_runs <- mse_runs$model_runs
extra_columns <- mse_runs$extra_columns2

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 130)

hcr_colors <- set_hcr_colors2(publication_hcrs)


### Spawning Biomass and Catch Plots
ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=publication_hcrs, om_filter=publication_oms)
plot_ssb(ssb_data, v1="hcr", v2="om", v3=NA, common_trajectory=common_trajectory, show_est = FALSE)
ggsave(filename=file.path(figures_dir, paste0("ssb", filetype)), width=width_small, height=height_small, units=c("in"))
 
plot_relative_ssb(ssb_data, v1="hcr", v2="om", common_trajectory = common_trajectory, base_hcr = "No Fishing")
ggsave(filename=file.path(figures_dir, paste0("rel_ssb", filetype)), width=width_small, height=height_small, units=c("in"))

catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)
plot_landed_catch(catch_data, v1="hcr", v2="om", common_trajectory = common_trajectory)
ggsave(filename=file.path(figures_dir, paste0("catch", filetype)), width=width_small, height=height_small, units=c("in"))

plot_ssb_catch(
    ssb_data,
    catch_data,
    v1="hcr",
    v2="om",
    common_trajectory = common_trajectory
)
ggsave(filename=file.path(here::here(), "figures", "ssb_catch.jpeg"), width=16, height=10, units="in")


### Performance Metrics
performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=publication_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb")
)
perf_data <- performance_metrics$perf_data %>% 
    mutate(
        name = factor(name, levels=publication_metrics)
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

plot_performance_metric_summary(perf_data, v1="hcr", v2="om")+
    theme(
        panel.spacing.x = unit(0.75, "cm"),
    )
ggsave(filename=file.path(here::here(), "figures", paste0("performance2", filetype)), width=18, height=12, units="in")

### Resilience Metrics
resilience_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=c("Immediate Crash Recruitment"),
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("crash_time", "recovery_time")
)

resilience_data <- resilience_metrics$perf_data %>%
    add_row(
        om = "Immediate Crash Recruitment",
        hcr = "No Fishing",
        .width = c(0.5, 0.8),
        .point = "median",
        .interval = "qi",
        name = "Crash Time",
        median = Inf, lower = Inf, upper = Inf
    ) %>%
    mutate(
        om = factor(om),
        hcr = factor(hcr, levels=publication_hcrs),
        name = factor(name, levels=c("Crash Time", "Recovery Time"))
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

plot_performance_metric_summary(resilience_data)+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("Year", limits=c(0, 20)),
            scale_x_continuous("Year", limits=c(0, 50))
        )
    )+
    theme(
        strip.text.y = element_blank(),
        plot.margin = margin(10, 30, 10, 10),
        panel.spacing.x = unit(30, "pt")
    )
ggsave(filename=file.path(here::here(), "figures", paste0("resilience2", filetype)), width=width_small, height=height_small, units="in")


### Aggregate Performance
full_performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=publication_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = FALSE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb") 
)

om_aggregated_performance <- full_performance_metrics$perf_data %>%
    group_by(sim, hcr) %>%
    summarise(
        across(2:6, \(x) mean(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(3:7, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_years"), 
            labels=publication_metrics
        )
    )

om_means <- om_aggregated_performance %>% group_by(name) %>% summarise(v = mean(value))

ggplot(om_aggregated_performance)+
    ggridges::stat_binline(aes(x=value, y=hcr, fill=hcr))+
    geom_vline(data=om_means, aes(xintercept=v), linetype="dashed")+
    scale_fill_manual(values=hcr_colors)+
    facet_wrap(~name, scales="free_x")+
    labs(y="")+
    custom_theme+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("", limits=c(0, 35)),
            scale_x_continuous("", limits=c(0, 0.06)),
            scale_x_continuous("", limits=c(0, 300)),
            scale_x_continuous("", limits=c(5, 10)),
            scale_x_continuous("", limits=c(0, 1), breaks=c(0, 0.25, 0.50, 0.75, 1.0), labels=seq(0, 1, 0.25), expand=c(0, NA))
        )
    )+
    guides(fill=guide_legend("Management \n Strategy", nrow=2))+
    theme(
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(1, "cm")
    )
ggsave(file.path(figures_dir, paste0("performance_metrics_agg_hist2", filetype)))


### Extra Analyses

# Get average performance metric value by HCR across the OM scenarios
om_aggregated_performance %>% group_by(hcr, name) %>%
    summarise(
        value = mean(value)
    ) %>%
    pivot_wider(names_from="name", values_from="value")

# Get number of simulations in which the population went extinct prematurely
ssb_data %>% filter(biomass < 1e-2) %>% pull(sim) %>% unique

# Get number of years in which stability and max harvest constraints were
# applied to their respective HCRs
catch_data %>%
    filter_times(time_horizon) %>%
    filter(fleet == "Fixed") %>%
    filter_hcr_om(hcrs = c("F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "25k Harvest Cap"), oms=publication_oms) %>%
    group_by(sim, om, hcr) %>%
    mutate(
        rel_perc_change = abs(total_catch/lag(total_catch)-1),
        constrained = case_when(
            hcr == "F40 +/- 5%" ~ ifelse(rel_perc_change >= 0.0499, 1, 0),
            hcr == "F40 +/- 10%" ~ ifelse(rel_perc_change >= 0.0999, 1, 0),
            hcr == "15k Harvest Cap" ~ ifelse(total_catch >= 14.99, 1, 0),
            hcr == "25k Harvest Cap" ~ ifelse(total_catch >= 24.99, 1, 0)
        )
    ) %>%
    # arrange(sim, om, time) %>%
    group_by(sim, om, hcr) %>%
    summarise(
        years_constrained = sum(constrained, na.rm=TRUE)/n()
    ) %>%
    group_by(om, hcr) %>%
    median_qi(years_constrained, .width = interval_widths) %>%
    print(n=100)

