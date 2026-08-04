rm(list=ls())
    
library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SPoRC)
library(tictoc)
library(doParallel)
library(afscOM)
library(patchwork)
library(afscOM) # may work but not certain

# Change to wherever your local copy of afscOM is
library(devtools)
# afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

# devtools::load_all(afscOM_dir)
# devtools::load_all(sablefishMSE_dir)

lapply(list.files("R", full.names = TRUE), source)

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons
nyears <- 110

sable_om <- readRDS("data/spatial_sablefish_om.RDS") # Read this saved OM from a file
# sable_om$model_options$fleet_apportionment <- matrix(c(0.80, 0.20), nrow=nrow(sable_om$model_options$fleet_apportionment), ncol=2, byrow=TRUE)

source("dev/oms.R")
source("dev/hcrs.R")


set.seed(895)
nsims <- 20
seed_list <- sample(1:(10*nsims), nsims)  # Draw 10 random seeds

mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 54
mse_options$recruitment_start_year <- 54
mse_options$n_proj_years <- 50
mse_options$run_estimation <- FALSE

mse_options_list <- listN(mse_options)

om_list <- listN(om_agemove_bhrec, om_agemove_regimerec, om_agemove_crashrec, om_climatemove_bhrec, om_climatemove_regimerec, om_climatemove_crashrec)
hcr_list <- listN(mp_f40_fullutil, mp_f50_fullutil, mp_f40hybrid_fullutil, mp_20cap_fullutil)#, mp_dynamicB0)#, mp_f40_fullutil, mp_f50, mp_f50_fullutil, mp_b30f50, mp_b30f50_fullutil, mp_f40hybrid, mp_f40hybrid_fullutil)

om_list <- listN(om_agemove_bhrec, om_agemove_regimerec, om_agemove_crashrec)
hcr_list <- listN(mp_f00_fullutil)

# mp5_modelrun <- run_mse_parallel(nsims, seed_list, om1, mp5, mse_options, nyears)
tic()
model_runs <- run_mse_multiple(
    om_list, 
    hcr_list, 
    seed_list,
    mse_options_list=mse_options_list,
    diagnostics = TRUE,
    save = TRUE
)
toc()
# extra_columns <- expand.grid(hcr=lapply(hcr_list, function(x) x$name), om=lapply(om_list, function(x) x$name))


# Data Processing
filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures", "publication")
width_small <- 12
height_small <- 8

width_large <- 16
height_large <- 16

width_medium <- 12
height_medium <- 12

# base_util_hcrs <- c("F40", "F50", "F50/B30", "F40 Hybrid")
# full_util_hcrs <- paste0(base_util_hcrs, " | Full Utilization")

publication_hcrs <- c("F40 | Full Utilization", "F50 | Full Utilization", "F40 Hybrid | Full Utilization", "20k Harvest Cap | Full Utilization")
publication_oms <- c("BH Recruit | AB Move", "BH Recruit | Climate Move", "Regime Recruit | AB Move", "Regime Recruit | Climate Move", "Crash Recruit | AB Move", "Crash Recruit | Climate Move")
publication_metrics <- c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years with Low SSB")

recruitment_scenarios <- c("BH Recruit", "Regime Recruit", "Crash Recruit")

extra_columns = expand.grid(
    om=publication_oms,
    hcr=publication_hcrs
)

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 105)

hcr_colors <- c("black", "#991c1c", "#001180", "#256c15", "#BB6A00", "#7500B0", "#5CABA3", "#BD4E98")
#names(hcr_colors) <- publication_hcrs

hcr_colors = c("black", "#256c15", "#0763d3", "#BD4E98", "#BB6A00")

hcr_filter <- publication_hcrs
om_filter <- publication_oms

###### Plot Diagnostics
convergence <- get_convergence_data(
    model_runs=NULL,
    extra_columns,
    hcr_filter,
    om_filter
)
write_csv(convergence, file.path("data", "zahneretal_2026_spatialmse_convergence.csv"))


unconverged <- convergence %>% filter(!converged)
converged_sims <- convergence %>% filter(converged) %>% pull(sim) %>% unique()

##### Load data #####
# Data files
ssb_data_file           <- file.path("data", "zahneretal_2026_spatialmse_SSB_biomass.csv")
catch_data_file         <- file.path("data", "zahneretal_2026_spatialmse_catch.csv")
econ_data_file          <- file.path("data", "zahneretal_2026_spatialmse_econvalue.csv")
avgage_catch_data_file  <- file.path("data", "zahneretal_2026_spatialmse_avgage_catch.csv")
avgage_pop_data_file    <- file.path("data", "zahneretal_2026_spatialmse_avgage_pop.csv")
napg_catch_data_file    <- file.path("data", "zahneretal_2026_spatialmse_napg_cat.csv")
napg_pop_data_file    <- file.path("data", "zahneretal_2026_spatialmse_napg_pop.csv")
abc_data_file           <- file.path("data", "zahneretal_2026_spatialmse_abctac.csv")
numbers_data_file       <- file.path("data", "zahneretal_2026_spatialmse_landed_numbers.csv")

ssb_data <- read_csv(ssb_data_file) %>% set_factor_levels
catch_data <- read_csv(catch_data_file) %>% set_factor_levels
econ_data <- read_csv(econ_data_file) %>% set_factor_levels
napg_catch_data <- read_csv(napg_catch_data_file) %>% set_factor_levels
napg_pop_data <- read_csv(napg_pop_data_file) %>% set_factor_levels
avgage_catch_data <- read_csv(avgage_catch_data_file) %>% set_factor_levels
avgage_pop_data <- read_csv(avgage_pop_data_file) %>% set_factor_levels
numbers_data <- read_csv(numbers_data_file) %>% set_factor_levels

#############################################
#############################################
#############################################
#### Spawning Biomass + Depletion ###########
#############################################
#############################################
#############################################
# ssb_data <- get_ssb_biomass(model_runs=NULL, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
# ssb_data <- ssb_data %>% separate_om_name %>% set_factor_levels
# write_csv(ssb_data, file.path("data", "zahneretal_2026_spatialmse_SSB_biomass.csv"))

make_ssb_plot <- function(ssb_data){
    reg_ssb_plots <- plot_ssb(
        ssb_data, 
        v1="hcr", v2="region", v3="om", show_est = FALSE,
        common_trajectory=common_trajectory 
    )+
    labs(title="Regional Spawning Biomass", y="SSB (1000s mt)")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 50), breaks=seq(0, 50, 10)),
            scale_y_continuous(limits=c(0, 50), breaks=seq(0, 50, 10)),
            scale_y_continuous(limits=c(0, 30), breaks=seq(0, 30, 10)),
            scale_y_continuous(limits=c(0, 125), breaks=seq(0, 125, 25)),
            scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 25))
        )
    )+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

    ssb_agg_plots <- plot_ssb(
            ssb_data, 
            v1="hcr", v2="om", v3=NA, show_est = FALSE,
            common_trajectory=common_trajectory
        )+
        coord_cartesian(ylim=c(0, 300), expand=0)+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )

    # ssb_agg_plots + depletion_plots + 
    reg_ssb_plots / ssb_agg_plots +
        plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
        theme(legend.position = "bottom")
}

make_ssb_plot(
    ssb_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter)
)
ggsave(filename=file.path(figures_dir, paste0("ssb_depletion", filetype)), width=width_large, height=height_large, units=c("in"))

ssb_data %>% separate_om_name() %>% set_factor_levels() %>%
plot_depletion(v1="hcr", v2="om", v3="region") +
    geom_hline(yintercept=c(0.35, 0.40), linetype="dashed", color="red")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 3.5), breaks=seq(0, 3.5, 0.5)),
            scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
            scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
            scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25)),
            scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))
        )
    )
ggsave(file.path(figures_dir, paste0("depletion", filetype)), width=width_large, height=height_large, units=c("in"))

ssb_data %>% filter_times(time_horizon) %>%
    filter(region != "Alaska") %>%
    mutate(
        Movement = factor(
            Movement, 
            levels=c("Base Move", "Climate Move"), 
            labels=c("Base", "Climate")
        )
    ) %>%
    group_by(Recruitment, Movement, hcr, region) %>%
    median_qi(spbio) %>%
    select(Recruitment, Movement, hcr, region, spbio) %>%

    ggplot()+
        geom_col(aes(y=spbio, x=Movement, fill=region))+
        scale_color_manual(values=hcr_colors)+
        facet_grid(Recruitment~hcr)+
        labs(title="Spawning Biomass by Movement Scenario", y="SSB (1000s mt)")+
        custom_theme

get_estimation_bias(ssb_data, om_em_cols=c("naa", "naa_est"), value_col="spbio", time_horizon) %>%
    plot_estimation_bias(type="boxplot")

#############################################
#############################################
#############################################
#### Landed Catch ###########################
#############################################
#############################################
#############################################
# catch_data <- get_landed_catch(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
# catch_data <- catch_data %>% separate_om_name %>% set_factor_levels
# write_csv(catch_data, file.path("data", "zahneretal_2026_spatialmse_catch.csv"))

# landed_num_data <- get_numbers(model_runs=NULL, extra_columns, dem_params=sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
# landed_num_data %>% separate_om_name %>% set_factor_levels %>% write_csv(file.path("data", "zahneretal_2026_spatialmse_landed_numbers.csv"))


make_catch_plot <- function(catch_data){

    reg_catch_plots <- plot_landed_catch(catch_data, v1="hcr", v2="om", v3="region", by_fleet=FALSE, common_trajectory = common_trajectory)+
        labs(title="Regional Landings", y="Catch (1000s mt)")+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2)),
                scale_y_continuous(limits=c(0, 8), breaks=seq(0, 8, 2)),
                scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 2)),
                scale_y_continuous(limits=c(0, 25), breaks=seq(0, 25, 5)),
                scale_y_continuous(limits=c(0, 15), breaks=seq(0, 15, 5))
            )
        )+
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )

    catch_agg_plot <- plot_landed_catch(catch_data, v1="hcr", v2="om", by_fleet=FALSE, common_trajectory = common_trajectory)+
        scale_y_continuous(limits=c(0, 60), name="Catch (1000s mt)")+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )

    reg_catch_plots / catch_agg_plot +  
        plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
        theme(legend.position = "bottom")
}

make_catch_plot(
    catch_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter)
)
ggsave(filename=file.path(figures_dir, paste0("catch", filetype)), width=width_large, height=height_large, units=c("in"))

catch_data %>% filter_times(time_horizon) %>%
    mutate(Movement = factor(Movement, levels=c("Base Move", "Climate Move"), labels=c("Base", "Climate"))) %>%
    group_by(Recruitment, Movement, hcr, region) %>%
    median_qi(total_catch) %>%
    select(Recruitment, Movement, hcr, region, total_catch) %>%

    ggplot()+
        geom_col(aes(y=total_catch, x=Movement, fill=region))+
        scale_color_manual(values=hcr_colors)+
        facet_grid(Recruitment~hcr)+
        labs(title="Total Catch by Movement Scenario", y="Total Catch (mt)")+
        custom_theme

#############################################
#############################################
#############################################
#### Dynamic Economic Value #################
#############################################
#############################################
#############################################
# econ_data <- get_dynamic_economic_value(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
# econ_data <- econ_data %>% separate_om_name %>% set_factor_levels
# write_csv(econ_data, file.path("data", "zahneretal_2026_spatialmse_econvalue.csv"))

make_econ_plot <- function(econ_data){

    reg_econ_plots <- plot_econ_value(econ_data, v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)+
        labs(title="Regional Economic Value", y="Economic Value")+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 30), breaks=seq(0, 30, 10)),
                scale_y_continuous(limits=c(0, 40), breaks=seq(0, 40, 10)),
                scale_y_continuous(limits=c(0, 30), breaks=seq(0, 30, 10)),
                scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 30)),
                scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 30))
                # scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 20)),
                # scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 20))
            )
        )+
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )

    econ_agg_plot <- plot_econ_value(econ_data, v1="hcr", v2="om", common_trajectory = common_trajectory)+
        scale_y_continuous(limits=c(0, 200))+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )

    reg_econ_plots / econ_agg_plot +  
        plot_layout(ncol=1, guides="collect", heights=c(1, 0.5))+
        plot_annotation(caption="Economic value currently calculated for Fixed gear fleet only") & 
        theme(legend.position = "bottom")

}

make_econ_plot(
    econ_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter)
)
ggsave(filename=file.path(figures_dir, paste0("economic_value", filetype)), width=width_large, height=height_large, units=c("in")) 

econ_data %>% filter_times(time_horizon) %>%
    mutate(Movement = factor(Movement, levels=c("Base Move", "Climate Move"), labels=c("Base", "Climate"))) %>%
    group_by(Recruitment, Movement, hcr, region) %>%
    median_qi(total_value) %>%
    select(Recruitment, Movement, hcr, region, total_value) %>%

    ggplot()+
        geom_col(aes(y=total_value, x=Movement, fill=region))+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(expand=0)+
        facet_grid(Recruitment~hcr)+
        labs(title="Total Economic Value by Movement Scenario", y="Total Value")+
        custom_theme
ggsave(filename=file.path(figures_dir, paste0("econ_region_props_movement", filetype)), width=width_large, height=height_small, units=c("in"))

#############################################
#############################################
#############################################
#### Recruitment ############################
#############################################
#############################################
#############################################
recruit_data <- get_recruits(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
# recruit_data <- recruit_data %>% separate_om_name %>% set_factor_levels %>%
#     mutate(region = as.factor(ifelse(is.na(region), "Alaska", as.character(region))))
# write_csv(recruit_data, file.path("data", "zahneretal_2026_spatialmse_recruitment.csv"))        

reg_rec_plots <- plot_recruitment(recruit_data, v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)
recruit_agg_plot <- plot_recruitment(recruit_data, v1="hcr", v2="om", show_est=FALSE, common_trajectory = common_trajectory)+
    scale_y_continuous("Recruits (millions)", limits=c(0, 100))+
    facet_wrap(~om, ncol=1)+
    labs(title="Alaska-Wide Recruitment")+
    theme(
        strip.background = element_blank(),
        strip.text = element_blank()
    )
recruit_agg_plot + reg_rec_plots + 
    plot_layout(nrow=1, guides="collect", axes = "collect", widths=c(0.5, 1)) & 
    theme(legend.position = "bottom")
ggsave(filename=file.path(figures_dir, paste0("recruitment", filetype)), width=16, height=8, units=c("in"))

get_estimation_bias(recruit_data, om_em_cols=c("naa", "naa_est"), value_col="rec", time_horizon) %>%
    plot_estimation_bias(type="boxplot")

get_estimation_bias(recruit_data, om_em_cols=c("naa", "naa_est"), value_col="rec", time_horizon) %>%
    group_by(om, hcr) %>%
    reframe(
        quants = quantile(bias, probs=c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0))
    ) %>%
    mutate(
        quantile = rep(c("min", "5%", "25%", "50%", "75%", "95%", "max"), length.out=7*3*2)
    ) %>%
    pivot_wider(names_from=quantile, values_from=quants)

recruit_data %>%
    filter_times(time_horizon) %>%
    filter(sim %in% sample(seed_list, 5)) %>%
    group_by(time, sim, om, Recruitment, Movement, hcr, L1) %>%
    summarise(rec = sum(rec, na.rm=TRUE)) %>%
    ggplot(aes(x=time, y=rec, color=hcr, fill=hcr, group=sim, alpha=0.8))+
        geom_line()+
        facet_grid(hcr~om)+
        custom_theme

plot_random_recruitment_trajectories(recruit_data, seed_list, 5, time_horizon)





recdevs <- process_big_outputs(NULL, var=c("global_rec_devs"), extra_columns, hcr_filter, om_filter, process_func=as.tibble)

recdevs %>% filter_times(time_horizon)

#############################################
#############################################
#############################################
#### ABC, TAC, Expected Landings ############
#############################################
#############################################
#############################################
# abc_tac_land <- get_management_quantities(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
# abc_tac_land <- abc_tac_land %>% separate_om_name %>% set_factor_levels
# write_csv(abc_tac_land, file.path("data", "zahneretal_2026_spatialmse_abctac.csv"))

abc_tac_land <- abctac_data %>% filter(fleet == "Trawl", L1 != "tac")

regional_abctac_plots <- plot_abc_tac(abc_tac_land, v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)+
    scale_color_manual(values=hcr_colors)+
    labs(title="Regional ABC and Expected Landings: Trawl Gear Fleet", y="ABC/Landings (1000s mt)")+
    ggh4x::facet_nested(
        rows=vars(region, L1), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 7)), scale_y_continuous(limits=c(0, 7)),
            scale_y_continuous(limits=c(0, 3)), scale_y_continuous(limits=c(0, 3)),
            scale_y_continuous(limits=c(0, 1)), scale_y_continuous(limits=c(0, 1)),
            scale_y_continuous(limits=c(0, 4)), scale_y_continuous(limits=c(0, 4)),
            scale_y_continuous(limits=c(0, 1)), scale_y_continuous(limits=c(0, 1))
        )
    )

agg_abctac_plot <- plot_abc_tac(abc_tac_land, v1="hcr", v2="om", v3=NA, common_trajectory = common_trajectory)+
    ggh4x::facet_nested(
        rows=vars(region, L1), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    guides(color="none", shape="none")+
    labs(y="ABC/Landings (1000s mt)")+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 50)),
            scale_y_continuous(limits=c(0, 50)),
            scale_y_continuous(limits=c(0, 50))
        )
    ) +
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

regional_abctac_plots / agg_abctac_plot +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5))+
    plot_annotation(caption="Expected landings accounts for ABC allocated to each region via apportionment and regional TAC utilization. Note: ABC = TAC != Expected Landings.") & 
    theme(legend.position = "bottom")
ggsave(file.path(here::here(), "figures", "abc_tac_land_trawl.jpeg"), width=width_large, height=height_large, units="in")



##### Average Age #####
# avgage_pop_data <- get_average_age(model_runs = NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter, pop_or_catch = "pop")
# avgage_pop_data <- avgage_pop_data %>% separate_om_name %>% set_factor_levels
# write_csv(avgage_pop_data, file.path("data", "zahneretal_2026_spatialmse_avgage_pop.csv"))

# avgage_catch_data <- get_average_age(model_runs = NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter, pop_or_catch = "catch")
# avgage_catch_data <- avgage_catch_data %>% separate_om_name %>% set_factor_levels
# write_csv(avgage_catch_data, file.path("data", "zahneretal_2026_spatialmse_avgage_catch.csv"))

make_avgage_plot <- function(average_age_data){

    reg_avgage_plots <- plot_average_age(average_age_data, v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)+
        labs(title="Average Population Age", y="Age (Years)")+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        scale_y_continuous(limits=c(0, 15))+
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )

    agg_avgage_plot <- plot_average_age(average_age_data, v1="hcr", v2="om", v3=NA, common_trajectory = common_trajectory)+
        # scale_color_manual(values=c("red", "blue", "green", "purple"))+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        guides(color="none", shape="none")+
        labs(y="Age (Years)")+
        scale_y_continuous(limits=c(0, 15))+
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )

    reg_avgage_plots / agg_avgage_plot +
        plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
        theme(legend.position = "bottom")
}


make_avgage_plot(
    avgage_pop_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter)
)
ggsave(file.path(here::here(), "figures", "average_age.jpeg"), width=width_large, height=height_large, units="in")


#############################################
#############################################
#############################################
#### Numbers at Price Grade #################
#############################################
#############################################
#############################################
# napg_pop_data <- get_numbers_at_age(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter, pop_or_catch = "pop")
# napg_pop_data <- napg_pop_data %>% separate_om_name %>% set_factor_levels
# write_csv(napg_pop_data, file.path("data", "zahneretal_2026_spatialmse_napg_pop.csv"))

# napg_catch_data <- get_numbers_at_age(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter, pop_or_catch = "catch")
# napg_catch_data <- napg_catch_data %>% separate_om_name %>% set_factor_levels
# write_csv(napg_catch_data, file.path("data", "zahneretal_2026_spatialmse_napg_catch.csv"))

make_napg_plot <- function(napg_data){
    reg_highpg_plots <- plot_numbers_at_pricegrade(napg_data, v1="hcr", v2="om", v3="region", price_grade="Grade 7+ (15+yo)", common_trajectory = common_trajectory)+
        labs(title="Average Number Grade 7+ Individuals", y="Individuals")+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        # scale_y_continuous(limits=c(0, 25))+
        ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 2.5)),
                scale_y_continuous(limits=c(0, 2.5)),
                scale_y_continuous(limits=c(0, 1.5)),
                scale_y_continuous(limits=c(0, 7.5)),
                scale_y_continuous(limits=c(0, 7.5))
            )
        )+
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )

    agg_highpg_plot <- plot_numbers_at_pricegrade(napg_data, v1="hcr", v2="om", v3=NA, price_grade="Grade 7+ (15+yo)", common_trajectory = common_trajectory)+
        ggh4x::facet_nested(
            rows=vars(region), 
            cols=vars(Recruitment, Movement), 
            scales="free_y"
        )+
        guides(color="none", shape="none")+
        labs(y="Number Grade 7+ Individuals")+
        scale_y_continuous(limits=c(0, 18))+
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )

    reg_highpg_plots / agg_highpg_plot +
        plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
        theme(legend.position = "bottom")
}

make_napg_plot(
    napg_catch_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter)
)
ggsave(file.path(here::here(), "figures", "napg7_catch.jpeg"), width=width_large, height=height_large, units="in")

make_napg_plot(
    napg_pop_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter)
)
ggsave(file.path(here::here(), "figures", "napg7_pop.jpeg"), width=width_large, height=height_large, units="in")


#################################
#################################
##### Fleet Proportions #########
#################################
#################################
fixedprop_catch_data <- catch_data %>%
    group_by(time, fleet, sim, om, hcr) %>% 
    summarise(
        alaska_catch = sum(total_catch),
        fleet_catch = sum(catch)
    ) %>%
    separate_om_name %>%
    filter(fleet=="Fixed") %>%
    mutate(prop = fleet_catch/alaska_catch)

fixedprop_econ_data <- econ_data %>%
    group_by(time, sim, region, om, hcr) %>% 
    mutate(
        fleet_value = sum(total_value),
    ) %>%
    group_by(time, fleet, sim, om, hcr) %>%
    summarise(
        alaska_value = sum(fleet_value),
        fleet_value = sum(total_value)
    ) %>%
    separate_om_name %>%
    filter(fleet=="Fixed") %>%
    mutate(prop = fleet_value/alaska_value)

# fixedprop_napg_data <- napg_catch_data %>%
#     filter(class == "Grade 7+ (15+yo)") %>%
#     group_by(time, sim, region, om, hcr) %>% 
#     mutate(
#         fleet_value = sum(value),
#     ) %>%
#     group_by(time, fleet, sim, om, hcr) %>%
#     summarise(
#         alaska_value = sum(fleet_value),
#         fleet_value = sum(value)
#     ) %>%
#     separate_om_name %>%
#     filter(fleet=="Fixed") %>%
#     mutate(prop = fleet_value/alaska_value)

fleet_props <- bind_rows(
    fixedprop_catch_data %>% mutate(metric="catch"),
    fixedprop_econ_data %>% mutate(metric="economic")
    # fixedprop_napg_data %>% mutate(metric="napg7")
) %>% select(time, fleet, sim, om, Recruitment, Movement, hcr, prop, metric) %>%
    group_by(time, om, Recruitment, Movement, hcr, fleet, metric) %>%
    median_qi(prop, .width = interval_widths)

ggplot(fleet_props %>% filter_times(time_horizon)) +
    geom_line(aes(x=time, y=prop, group=interaction(fleet, hcr), color=hcr), size=0.85)+
    scale_color_manual(values=hcr_colors)+
    ggh4x::facet_nested(
        cols=vars(metric), 
        rows=vars(Recruitment, Movement)
    )+
    labs(title="Proportion of Catch and Economic Value for Fixed Gear Fleet", y="Proportion of Catch/Economic Value", x="Year")+
    custom_theme

ggsave(file.path(figures_dir, paste0("prop_fixed_gear", filetype)), width=width_medium, height=height_medium, units=c("in"))

############################################
############################################
##### Proportions of Metrics by Region #####
############################################
############################################

ssb_props <- ssb_data %>% select(-biomass, -ssb0, -dep) %>% filter(L1 == "naa") %>%
    calculate_regional_props(val_name="spbio")
# catch_props <- catch_data %>% filter(fleet == "Fixed") %>% select(-fleet, -catch) %>%
#     calculate_regional_props(val_name="total_catch")
# econ_props <- econ_data %>% group_by(time, sim, region, om, hcr) %>% 
#     mutate(total_value = sum(total_value)) %>%
#     calculate_regional_props(val_name="total_value")
# napg_props <- napg_catch_data %>% filter(class == "Grade 7+ (15+yo)") %>%
#     calculate_regional_props(val_name="value")
props <- bind_rows(
    ssb_props %>% mutate(metric="ssb"),
    catch_props %>% mutate(metric="catch"),
    econ_props %>% mutate(metric="economic_value"),
    napg_props %>% mutate(metric="napg7")
) %>% 
    select(Recruitment, Movement, hcr, region, name, median, lower, upper, .width, metric) %>%
    # set_factor_levels %>%
    mutate(
        metric = factor(
            metric, 
            levels=c("ssb", "catch", "napg7", "economic_value"), 
            labels=c("SSB", "Catch", "High Price\nIndividuals", "Economic Value")
        ),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA", "Alaska"))
    )

props <- read_csv(file.path("data", "zahneretal_2026_spatialmse_regionalprops.csv")) %>% set_factor_levels

props_regional_plot <- props %>%
    filter(
        name == "prop", .width==0.50,
        Recruitment == "BH Recruit", 
        hcr %in% c("F40 | Full Utilization", "F40 Hybrid | Full Utilization")
    ) %>%
    ggplot()+
        geom_bar(aes(x=median, y=metric, fill=region), stat="identity", position="fill")+
        scale_color_manual(values=hcr_colors)+
        scale_x_continuous(labels=seq(0, 1, 0.25))+
        geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype="dashed", color="black")+
        coord_cartesian(expand=0)+
        facet_grid(Movement~hcr)+
        labs(title="Regional Proportions of Key Metrics", y="Movement Scenario", x="Proportion")+
        custom_theme+
        theme(
            panel.spacing.x = unit(35, "pt"),
        )
ggsave(file.path(figures_dir, paste0("regional_props_allmetrics", filetype)), width=width_large, height=height_medium, units=c("in"))

delta_props_regional_plot <- props %>%
    filter(
        name == "scaled_delta",
        Recruitment == "BH Recruit", 
        hcr %in% c("F40 | Full Utilization", "F40 Hybrid | Full Utilization")
    ) %>%
    # pivot_wider(names_from="Movement", values_from="prop") %>%
    # mutate(delta = (`Climate Move`-`Base Move`)/`Base Move`) %>%

    ggplot()+
        geom_bar(aes(y=median, x=metric, fill=region), stat="identity", position="dodge")+
        geom_pointinterval(aes(x=metric, y=median, ymin=lower, ymax=upper, color=region), point_size=0.01, position=position_dodge(width=0.9))+
        # scale_color_manual(values=hcr_colors)+
        # scale_y_continuous(limits=c(-0.25, 0.25), breaks=seq(-0.25, 0.25, 0.0625), labels=seq(-0.25, 0.25, 0.0625))+
        scale_y_continuous(limits=c(-0.75, 5.5), breaks=c(0, seq(-0.75, 5.5, 0.75)))+
        scale_x_discrete(labels=c("SSB", "Catch", "High Price\nIndividuals", "Economic\nValue"))+
        # scale_y_continuous(limits=c(-0.25, 0.25), breaks=seq(-0.25, 0.25, 0.125), labels=seq(-0.25, 0.25, 0.125))+
        # geom_vline(xintercept = c(0.25, 0.5, 0.75), linetype="dashed", color="black")+
        # coord_cartesian(expand=0, ylim=c(-0.3, 0.25))+
        guides(color="none")+
        facet_wrap(~hcr)+
        labs(title="Change in Regional Proportions of Key Metrics", y="Change in Proportion\n(Base Move -> Climate Move)", x="Metric")+
        custom_theme+
        theme(
            panel.spacing.x = unit(35, "pt"),
        )


props_regional_plot / delta_props_regional_plot + 
    plot_layout(ncol=1, guides="collect") & theme(legend.position="bottom")
ggsave(file.path(figures_dir, paste0("regional_props_allmetrics_scaled", filetype)), width=width_large, height=height_large, units="in")


##############################################################
##############################################################
##### Decomposition of SSB and Catch by Age and Biomass ######
##############################################################
##############################################################

ssb_data <- read_csv(ssb_data_file) %>% set_factor_levels
avgage_pop_data <- read_csv(avgage_pop_data_file) %>% set_factor_levels
catch_data <- read_csv(catch_data_file) %>% set_factor_levels
avgage_catch_data <- read_csv(avgage_catch_data_file)  %>% set_factor_levels
avgprice_data <- get_average_dynamic_price(model_runs=NULL, extra_columns, hcr_filter, om_filter)

age_ssb_decomp_data <- numbers_data %>% 
    # filter(L1=="naa") %>% 
    # select(-L1) %>%
    calculate_decomp(avgage_pop_data, val_name="total_num") %>%
    mutate(
        time=factor(time),
        name = factor(name, levels=c("quality_effect", "quantity_effect", "net_effect"), labels=c("Age Effect", "Population Size Effect", "Net SSB Effect"))
    )

age_ssb_decomp_data %>% #filter(hcr == "F40 | Full Utilization") %>%
    filter(Recruitment == "BH Recruit") %>%
    ggplot()+
        geom_bar(aes(x=time, y=median, fill=hcr), width=0.75, position="dodge", stat="identity")+
        geom_pointrange(aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), size=0.01, position=position_dodge(width=0.75))+
        scale_x_discrete(breaks=seq(55, 105, 5))+
        ggh4x::facet_nested(
            rows=vars(name, Movement), 
            cols=vars(region),
            scales="free_y"
        )+
        labs(title="Spawning Biomass Decomposition under BH Recruit", y="Effect Size", x="Year")+
        custom_theme+
        theme(
            # panel.spacing.x = unit(rep_len(c(0.1, 0.5), 7), "in")
            panel.spacing.y = unit(rep_len(c(0.1, 0.5), 5), "in")
        )
ggsave(file.path(figures_dir, paste0("biomass_decomp_bhrecruit_hcr", filetype)), width=width_large, height=height_large, units="in")


age_catch_decomp_data <- numbers_data %>% 
    # separate_om_name %>% set_factor_levels %>%
    # filter(L1=="land_caa") %>% 
    # select(-L1) %>%
    # group_by(time, sim, om, Recruitment, Movement, hcr, region) %>%
    # summarise(total_catch = sum(catch, na.rm=TRUE)) %>%
    calculate_decomp(avgage_catch_data, val_name="total_dead") %>%
    mutate(
        time=factor(time),
        name = factor(name, levels=c("quality_effect", "quantity_effect", "net_effect"), labels=c("Quality Effect", "Quantity Effect", "Net Effect"))
    )

age_catch_decomp_data %>% #filter(hcr == "F40 | Full Utilization") %>%
    filter(Recruitment == "BH Recruit") %>%
    ggplot()+
        geom_bar(aes(x=time, y=median, fill=hcr), width=0.75, position="dodge", stat="identity")+
        geom_pointrange(aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), size=0.01, position=position_dodge(width=0.75))+
        scale_x_discrete(breaks=seq(55, 105, 5))+
        ggh4x::facet_nested(
            rows=vars(name, Movement), 
            cols=vars(region),
            scales="free_y"
        )+
        labs(title="Catch Decomposition under BH Recruitment Scenario", y="Effect Size", x="Year")+
        custom_theme+
        theme(
            panel.grid.major.x = element_blank(),
            panel.spacing.y = unit(rep_len(c(0.1, 0.5), 5), "in")
        )
ggsave(file.path(figures_dir, paste0("catch_decomp_bhercruit_hcr", filetype)), width=width_large, height=height_large, units="in")

price_catch_decomp_data <- catch_data %>% 
    filter(L1=="land_caa", fleet=="Fixed") %>% 
    select(-L1) %>%
    calculate_decomp(avgprice_data %>% separate_om_name %>% set_factor_levels %>% rename(avg_age=average_price), val_name="total_catch") %>%
    mutate(
        time=factor(time),
        name = factor(name, levels=c("quality_effect", "quantity_effect", "net_effect"), labels=c("Price Effect", "Volume Effect", "Net Effect"))
    )

price_catch_decomp_data %>% #filter(hcr == "F40 | Full Utilization") %>%
    filter(Recruitment == "BH Recruit") %>%
    ggplot()+
        geom_bar(aes(x=time, y=median, fill=hcr), width=0.75, position="dodge", stat="identity")+
        geom_pointrange(aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), size=0.01, position=position_dodge(width=0.75))+
        scale_x_discrete(breaks=seq(55, 105, 5))+
        ggh4x::facet_nested(
            rows=vars(name, Movement), 
            cols=vars(region),
            scales="free_y"
        )+
        labs(title="Age-Biomass Decomposition under BH Recruitment Scenario", y="Effect Size", x="Year")+
        custom_theme+
        theme(
            # panel.spacing.x = unit(rep_len(c(0.1, 0.5), 7), "in")
            panel.spacing.y = unit(rep_len(c(0.1, 0.5), 5), "in")
        )
ggsave(file.path(figures_dir, paste0("value_decomp_bhrecruit_effect", filetype)), width=width_large, height=height_large, units="in")


#####################################
#####################################
### Comparison to Single Are MSE ####
#####################################
#####################################

agg_ssb_data <- ssb_data %>% separate_om_name %>% filter(L1=="naa", Movement=="AB Move") %>% select(-biomass, -Movement) %>%
    group_by(time, sim, L1, Recruitment, hcr) %>%
    summarise(
        spbio = sum(spbio, na.rm=TRUE)
    ) %>%
    group_by(time, Recruitment, hcr) %>%
    median_qi(spbio, .width = interval_widths) %>%
    separate(hcr, sep=c("\\s[|]\\s"), into=c("hcr", "utilization"), remove=FALSE) %>%
    mutate(
        model = "Spatial"
    ) %>%
    select(-utilization)

agg_catch_data <- catch_data %>% filter(fleet == "Fixed", Movement=="Base Move") %>% select(-Movement) %>%
    group_by(time, sim, L1, Recruitment, hcr) %>%
    summarise(
        total_catch = sum(total_catch, na.rm=TRUE)
    ) %>%
    group_by(time, Recruitment, hcr) %>%
    median_qi(total_catch, .width = interval_widths) %>%
    separate(hcr, sep=c("\\s[|]\\s"), into=c("hcr", "utilization"), remove=FALSE) %>%
    mutate(
        model = "Spatial"
    ) %>%
    select(-utilization)

singlearea_ssb_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_agg_SSB_biomass.csv")) %>%
    select(-biomass, -region) %>% filter(L1=="naa") %>% 
    rename(Recruitment=om) %>%
    group_by(time, Recruitment, hcr) %>%
    median_qi(spbio, .width = interval_widths) %>%
    mutate(
        model = "Single Area",
        Recruitment = factor(Recruitment, levels=c("Beverton-Holt Recruitment", "Regime Recruitment", "Crash Recruitment"), labels=c("BH Recruit", "Regime Recruit", "Crash Recruit"))
    )

singlearea_catch_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_agg_landed_catch.csv")) %>%
    filter(fleet=="Fixed") %>% 
    select(-catch) %>%
    rename(Recruitment=om) %>%
    group_by(time, Recruitment, hcr) %>%
    median_qi(total_catch, .width = interval_widths) %>%
    mutate(
        model = "Single Area",
        Recruitment = factor(Recruitment, levels=c("Beverton-Holt Recruitment", "Regime Recruitment", "Crash Recruitment"), labels=c("BH Recruit", "Regime Recruit", "Crash Recruit"))
    )

agg_recruit_data <- recruit_data %>% separate_om_name %>% filter(L1=="naa", Movement=="AB Move") %>% select(-Movement) %>%
    group_by(time, sim, L1, Recruitment, hcr) %>%
    summarise(
        recruits = sum(rec, na.rm=TRUE)
    ) %>%
    group_by(time, Recruitment, hcr) %>%
    median_qi(recruits, .width = interval_widths) %>%
    separate(hcr, sep=c("\\s[|]\\s"), into=c("hcr", "utilization"), remove=FALSE) %>%
    mutate(
        model = "Spatial"
    ) %>%
    select(-utilization)

singlearea_recruit_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_agg_recruits.csv")) %>%
    rename(Recruitment=om, recruits=rec) %>%
    group_by(time, Recruitment, hcr) %>%
    median_qi(recruits, .width = interval_widths) %>%
    mutate(
        model = "Single Area",
        Recruitment = factor(Recruitment, levels=c("Beverton-Holt Recruitment", "Regime Recruitment", "Crash Recruitment"), labels=c("BH Recruit", "Regime Recruit", "Crash Recruit"))
    )

agg_ssb_data %>%
    bind_rows(singlearea_ssb_data) %>%
    ggplot()+
        geom_line(aes(x=time, y=spbio, group=interaction(model, hcr), color=hcr, linetype=model), size=0.85)+
        scale_y_continuous(limits=c(0, 420), breaks=seq(0, 420, 50))+
        scale_color_manual(values=hcr_colors)+
        facet_wrap(~Recruitment)+
        labs(title="Alaska-Wide Spawning Biomass", y="SSB (1000s mt)", x="Year")+
        custom_theme

agg_catch_data %>%
    bind_rows(singlearea_catch_data) %>%
    ggplot()+
        geom_line(aes(x=time, y=total_catch, group=interaction(model, hcr), color=hcr, linetype=model), size=0.85)+
        scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10))+
        scale_color_manual(values=hcr_colors)+
        facet_wrap(~Recruitment)+
        labs(title="Alaska-Wide Catch", y="Catch (1000s mt)", x="Year")+
        custom_theme

agg_recruit_data %>%
    bind_rows(singlearea_recruit_data) %>%
    ggplot()+
        geom_line(aes(x=time, y=recruits, group=interaction(model, hcr), color=hcr, linetype=model), size=0.85)+
        scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 20))+
        scale_color_manual(values=hcr_colors)+
        facet_wrap(~Recruitment)+
        labs(title="Alaska-Wide Recruitment", y="Recruits (1000s)", x="Year")+
        custom_theme

#####################################
#####################################
###### Performance Metrics ##########
#####################################
#####################################

publication_metrics <- c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion Large Catch", "Dynamic Annual Value")
publication_metric_names <- c("Annual\nCatch", "Catch\nAAV", "SSB", "Average\nAge", "Average\nCatch\nLarge", "Dynamic\nValue")

tic()
performance_metrics_regional <- performance_metric_summary(
    model_runs=NULL, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr", "region"),
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "avg_catch_lg", "dynamic_value")
)

performance_metrics_alaska <- performance_metric_summary(
    model_runs=NULL, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "avg_catch_lg", "dynamic_value")
)
toc()

performance_metrics_regional <- readRDS(file.path("data", "performance_metrics_regional.RDS"))
performance_metrics_alaska <- readRDS(file.path("data", "performance_metrics_alaska.RDS"))

perf_data <- bind_rows(
    performance_metrics_regional$perf_data,
    performance_metrics_alaska$perf_data %>% mutate(region="Alaska")
) %>%
    filter_hcr_om(publication_hcrs, publication_oms) %>%
    filter(.width == 0.5, region != "Alaska") %>%
    group_by(om, region, name) %>%
    scale_and_rank("median") %>%
    mutate(
            hcr = factor(hcr, levels=hcr_filter),
            om = factor(om, levels=om_filter),
            name = factor(name, levels=publication_metrics, labels=publication_metric_names)
    ) %>%
    separate_om_name %>%
    set_factor_levels

# Performance Radar Plots
perf_data %>% 
    select(om, hcr, region, name, scaled, Recruitment, Movement) %>% 
    arrange(om, region, name) %>%
    mutate(
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA", "Alaska"))
    ) %>%
    # filter(Movement == "Base Move") %>%
    ggplot(aes(x=name, y=scaled, color=hcr, fill=hcr, group=hcr))+
        geom_point(size=3)+
        geom_line()+
        geom_polygon(alpha=0)+
        scale_y_continuous(breaks=seq(0, 1, 0.125))+
        # coord_cartesian(ylim=c(0.5, 1))+
        scale_color_manual(values=hcr_colors)+
        ggiraphExtra::coord_radar()+
        ggh4x::facet_nested(
            cols=vars(Recruitment, Movement), 
            rows=vars(region)
        )+
        custom_theme+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.spacing.x = grid::unit(rep_len(c(0.1, 0.5), 5), "in")
            # panel.spacing.y = unit(c(0.1, 0.1, 0.1, 0.1, 0.5), "in")        
        )
ggsave(file.path(figures_dir, paste0("performance_radar_omstd_hcr", filetype)), width=width_large, height=height_large, units="in")

performance_metrics_alaska_full <- performance_metric_summary(
    model_runs=NULL, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr", "region"),
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    summary_out = FALSE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "avg_catch_lg", "dynamic_value")
)
performance_metrics_alaska_full <- readRDS(file.path("data", "performance_metrics_alaska_full.RDS"))

perf_data_ak_full <- performance_metrics_alaska_full$perf_data %>% 
    mutate(region="Alaska") %>%
    filter_hcr_om(publication_hcrs, publication_oms) %>%
    separate_om_name %>%
    set_factor_levels %>%
    ungroup() %>%
    pivot_longer(6:11, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_large", "dyn_annual_value"), 
            labels=publication_metrics
        )
    ) %>%
    group_by(sim, om, Recruitment, Movement, region, name) %>%
    scale_and_rank("value")

# Calculate radar area
perf_data_ak_full %>% select(om, hcr, region, name, scaled) %>% arrange(om, region, name) %>%
    mutate(name = factor(name, levels=c(publication_metrics[1], rev(publication_metrics[2:6])), labels=publication_metrics)) %>%
    group_by(om, region) %>%
    group_split() %>%
    purrr::map(function(df){

        df %>% group_by(hcr) %>%
            mutate(
                weight = case_when(name == publication_metrics[1] ~ 1.5,
                                   name == publication_metrics[2] ~ 1.0,
                                   name == publication_metrics[3] ~ 1.5,
                                   name == publication_metrics[4] ~ 2.0,
                                   name == publication_metrics[5] ~ 2.5,
                                   name == publication_metrics[6] ~ 1.5),
                scaled = scaled * weight,
                n = row_number()-1,
                theta = 360/n(),
                n_theta = n * theta,
                theta_rad = theta * pi/180,
                n_theta_rad = n * theta_rad,
                x = scaled * cos(n_theta_rad),
                y = scaled * sin(n_theta_rad),
            ) %>%
            select(om, hcr, region, name, scaled, Recruitment, Movement, x, y) %>%
            group_by(hcr) %>%
            mutate(
                area = pracma::polyarea(as.numeric(x), as.numeric(y))
            ) %>%
            distinct(om, hcr, region, Recruitment, Movement, area)
        
    }) %>%
    bind_rows() %>%
    pivot_wider(names_from=hcr, values_from=area) %>%
    arrange(om) %>%
    print(n=100)

topsis_alaska_full <- compute_topsis(
    perf_data = perf_data_ak_full, #%>% rename(value=median),
    topsis_splits = c("sim", "om"),
    topsis_weights = c(0.2, 0.1, 0.15, 0.1, 0.2, 0.25),
    topsis_minmax = c("max", "min", "max", "max", "max", "max")
) %>%
    pivot_longer(3:6, names_to="hcr", values_to="topsis_score") %>%
    # pivot_wider(names_from="region", values_from="topsis_score") %>%
    separate_om_name %>%
    group_by(Recruitment, Movement, hcr) %>%
    median_qi(topsis_score, .width=interval_widths) %>%
    mutate(region="Alaska") %>%
    set_factor_levels

ggplot(topsis_alaska_full)+
    geom_pointinterval(aes(x=topsis_score, xmin=.lower, xmax=.upper, y=hcr, color=hcr))+
    scale_color_manual(values=hcr_colors)+
    facet_grid(Movement~Recruitment)+
    custom_theme
ggsave(file.path(figures_dir, paste0("topsis_ak", filetype)), width=width_small, height=height_medium, units="in")

performance_metrics_regional_full <- performance_metric_summary(
    model_runs=NULL, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr", "region"),
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    summary_out = FALSE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "avg_catch_lg", "dynamic_value")
)
performance_metrics_regional_full <- readRDS(file.path("data", "performance_metrics_regional_full.RDS"))

perf_data_reg_full <- performance_metrics_regional_full$perf_data %>% 
    filter_hcr_om(publication_hcrs, publication_oms) %>%
    separate_om_name %>%
    set_factor_levels %>%
    ungroup() %>%
    pivot_longer(7:12, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_large", "dyn_annual_value"), 
            labels=publication_metrics
        )
    ) %>%
    group_by(sim, om, Recruitment, Movement, region, name) %>%
    scale_and_rank("value")
    
topsis_regional_full <- compute_topsis(
    perf_data = perf_data_reg_full, #%>% rename(value=va),
    topsis_splits = c("sim", "om", "region"),
    topsis_weights = rep(1/6, 6), #c(0.2, 0.1, 0.15, 0.1, 0.2, 0.25),
    # topsis_weights = c(0.15, 0.1, 0.15, 0.1, 0.2, 0.3),
    topsis_minmax = c("max", "min", "max", "max", "max", "max")
) %>%
    pivot_longer(4:7, names_to="hcr", values_to="topsis_score") %>%
    # pivot_wider(names_from="region", values_from="topsis_score") %>%
    group_by(om, hcr, region) %>%
    median_qi(topsis_score, .width = interval_widths) %>%
    separate_om_name %>%
    set_factor_levels

    # print(n=100)
ggplot(topsis_regional_full)+
    geom_pointinterval(aes(x=topsis_score, xmin=.lower, xmax=.upper, y=hcr, color=hcr), position=position_dodge(width=0.5))+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement)
    )+
    scale_color_manual(values=hcr_colors)+
    labs(title="Topsis Score by HCR and Recruitment Scenario", y="HCR", x="Topsis Score")+
    custom_theme
ggsave(file.path(figures_dir, paste0("topsis_region", filetype)), width=width_large, height=height_large, units="in")

topsis %>%
    pivot_longer(4:7, names_to="hcr", values_to="topsis_score") %>%
    # pivot_wider(names_from="region", values_from="topsis_score") %>%
    group_by(om, hcr, region) %>%
    median_qi(topsis_score, .width = interval_widths) %>%
    separate_om_name %>%
    set_factor_levels





##### Resiliency Metrics #####
crash_time <- biomass_crash_time(
    model_runs, 
    extra_columns2, 
    sable_om$dem_params, 
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    hcr_filter=hcr_names,
    om_filter=c("Crash Recruitment")
)

recovery_time <- biomass_recovery_time(
    model_runs, 
    extra_columns2, 
    sable_om$dem_params, 
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    hcr_filter=hcr_names,
    om_filter=c("Crash Recruitment")
)


crash_recovery_time <- crash_time %>% reformat_ggdist_long(n=2) %>%
     bind_rows(recovery_time %>% reformat_ggdist_long(n=2)) %>%
     filter(hcr %in% publication_hcrs) %>%
     mutate(
            name=factor(name, labels=c("Crash Time", "Recovery Time")),
            hcr = factor(hcr, levels=publication_hcrs)
    )

ggplot(crash_recovery_time)+
    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=hcr, color=hcr), point_size=3, position="dodge")+
    facet_wrap(~name, scales="free_x")+
    # geom_vline(data=summary, aes(xintercept = median), color="black")+
    scale_shape_discrete()+
    scale_color_manual(values=hcr_colors)+
    scale_x_continuous("Years", limits=c(0, 20))+
    # facet_wrap(vars(name), scales="free_x")+
    labs(y="", x="", shape="OM", color="HCR")+
    coord_cartesian(expand=0)+
    guides(shape=guide_legend(nrow=3), color="none")+
    custom_theme+
    theme(
        plot.title = element_text(size=16),
        plot.margin = margin(0.24, 1, 0.25, 0.25, unit="cm"),
        panel.spacing.x = unit(1, "cm"),
    )
ggsave(file.path(here::here(), "figures", "resilience_metrics.jpeg"), width=12, height=8, units="in")

