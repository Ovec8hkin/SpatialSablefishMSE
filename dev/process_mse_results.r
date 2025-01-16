library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggdist)
library(broom)
library(afscOM)

lapply(list.files("R", full.names = TRUE), source)

filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures")

nyears <- 110
sable_om <- readRDS("data/sablefish_om_big.RDS")

# source(file.path(here::here(), "dev", "oms.R"))
# source(file.path(here::here(), "dev", "hcrs.R"))
# om_list <- afscOM::listN(om_rand_recruit, om_bh_recruit, om_bhcyclic_recruit, om_immcrash_recruit)
# hcr_list <- afscOM::listN(
#     mp_f40, mp_f50, mp_b30f40, mp_b40f50,
#     mp_5perc, mp_10perc, mp_10perc_up, mp_15perc,
#     mp_15cap, mp_25cap,
#     mp_f50chr,
#     mp_pfmc4010, mp_bcsable,
#     mp_f00chr
# )

# om_names <- unlist(lapply(om_list, \(x) x$name))
# hcr_names <- unlist(lapply(hcr_list, \(x) x$name))

publication_hcrs <- c("F40", "F50", "B40/F50", "F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "25k Harvest Cap", "Constant F50", "PFMC 40-10", "British Columbia", "No Fishing")
publication_oms <- c("Random Recruitment", "Beverton-Holt Cyclic Recruitment", "Immediate Crash Recruitment")
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years with Low SSB")

mse_runs <- get_saved_model_runs(om_order=publication_oms, hcr_order=publication_hcrs)
model_runs <- mse_runs$model_runs
extra_columns <- mse_runs$extra_columns2

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 130)

width_small <- 12
height_small <- 8

hcr_colors <- set_hcr_colors2(publication_hcrs)
rank_colors <- c(
    "#D55E00",
    "#FF740A",
    "#FF8B33",
    "#FFA35C",
    "#FFBA85",
    "#AAAAAA",
    "#5CC3FF",
    "#33B4FF",
    "#0AA5FF",
    "#008EE0",
    "#0072B2"
)


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
write_csv(performance_metrics$perf_data, file=file.path(here::here(), "data", "perfs.csv"))

perf_data <- performance_metrics$perf_data %>% 
    mutate(
        name = factor(name, levels=publication_metrics)
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

plot_performance_metric_summary(perf_data, v1="hcr", v2="om")
ggsave(filename=file.path(here::here(), "figures", paste0("performance", filetype)), width=18, height=12, units="in")
 
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
        strip.text.y = element_blank()
    )
ggsave(filename=file.path(here::here(), "figures", paste0("resilience", filetype)), width=width_small, height=height_small, units="in")

### Aggregate Performance Distributions

perf_tradeoffs <- performance_metric_summary(
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

tradeoff_summ <- perf_tradeoffs$perf_data %>%
    # filter(hcr != "15k Harvest Cap") %>%
    group_by(sim, hcr) %>%
    summarise(
        across(2:6, \(x) median(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(3:7, names_to="name", values_to="value") %>%
    mutate(
        name = factor(name, levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_years"), labels=publication_metrics)
    )
    # filter(!(hcr == "No Fishing" & name %in% c("Annual Catch", "Catch AAV")))

ggplot(tradeoff_summ, aes(x=value, y=hcr, fill=hcr))+
    stat_slab(normalize="groups") +
    geom_point(aes(y="15k Harvest Cap", x=15), color=hcr_colors[names(hcr_colors) == "15k Harvest Cap"])+
    geom_point(aes(y="No Fishing", x=0), color=hcr_colors[names(hcr_colors) == "No Fishing"])+
    scale_fill_manual(values=hcr_colors)+
    facet_wrap(~name, scales="free_x")+
    # coord_cartesian(expand=0)+
    custom_theme+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("", limits=c(-1, 35)),
            scale_x_continuous("", limits=c(-0.001, 0.06)),
            scale_x_continuous("", limits=c(50, 300)),
            scale_x_continuous("", limits=c(5, 10)),
            scale_x_continuous("", limits=c(0.01, 1), breaks=c(0.03, 0.25, 0.50, 0.75, 1.0), labels=seq(0, 1, 0.25), expand=c(0, NA))
        )
    )+
    theme(
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(1, "cm")
    )

tradeoff_summ %>% filter(hcr == "No Fishing", name == "Annual Catch")
