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
nsims <- 200
seed_list <- sample(1:(10*nsims), nsims)  # Draw 10 random seeds

mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 54
mse_options$recruitment_start_year <- 54
mse_options$n_proj_years <- 50
mse_options$run_estimation <- TRUE

mse_options_list <- listN(mse_options)

om_list <- listN(om_agemove_bhrec, om_agemove_regimerec, om_agemove_crashrec, om_climatemove_bhrec, om_climatemove_regimerec, om_climatemove_crashrec)
hcr_list <- listN(mp_f40_fullutil, mp_f50_fullutil, mp_f40hybrid_fullutil, mp_20cap_fullutil)#, mp_dynamicB0)#, mp_f40_fullutil, mp_f50, mp_f50_fullutil, mp_b30f50, mp_b30f50_fullutil, mp_f40hybrid, mp_f40hybrid_fullutil)

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
figures_dir <- file.path(here::here(), "figures")
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
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years with Low SSB")

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

# plot_diags(
#     hcr = "F40Test",
#     om = "Regime Recruit | AB Move",
#     simulation_seed = 1602,
#     simulation_year = 1,
#     om_list, 
#     seed_list
# )

#############################################
#############################################
#############################################
#### Spawning Biomass + Depletion ###########
#############################################
#############################################
#############################################
ssb_data <- get_ssb_biomass(model_runs=NULL, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)

ssb_data <- ssb_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )

write_csv(ssb_data, file.path("data", "zahneretal_2026_spatialmse_SSB_biomass.csv"))

get_estimation_bias(ssb_data, om_em_cols=c("naa", "naa_est"), value_col="spbio", time_horizon) %>%
    plot_estimation_bias(type="timeseries")

ssb_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_SSB_biomass.csv")) %>% mutate(
    Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
    Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
    region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
) #%>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == "BH Recruit")

for(r in recruitment_scenarios){
    fname <- paste0("ssb_depletion_fullutil_", sub(" ", "_", tolower(r)))
    make_ssb_plot(
        ssb_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

for(r in recruitment_scenarios){
    fname <- paste0("ssb_depletion_baseutil_", sub(" ", "_", tolower(r)))
    make_ssb_plot(
        ssb_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)# %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

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
        coord_cartesian(ylim=c(0, 300))+
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


ggsave(filename=file.path(figures_dir, paste0("ssb_depletion", filetype)), width=width_large, height=height_large, units=c("in"))


ssb_props <- ssb_data %>% select(-biomass) %>% filter(L1 == "naa") %>%
    pivot_wider(names_from=region, values_from=spbio) %>%
    rowwise() %>%
    mutate(alaska = sum(across(BS:EGOA))) %>%
    mutate(across(BS:EGOA, ~ .x/alaska)) %>%
    filter_times(time_horizon) %>%
    pivot_longer(BS:EGOA, names_to="region", values_to="prop") %>%
    group_by(om, hcr, region) %>%
    median_qi(prop)

ssb_props %>% select(om, hcr, region, prop) %>%
    mutate(
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA")),
        om = factor(om, levels=rev(publication_oms))
    ) %>%
    # pivot_wider(names_from=region, values_from=prop) %>%
    # arrange(hcr) %>%

    ggplot()+
        geom_bar(aes(y=om, x=prop, fill=region), stat="identity", position="fill")+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(expand=0)+
        facet_wrap(~hcr, ncol=2)+
        custom_theme

ssb_data %>% filter_times(time_horizon) %>%
    mutate(Movement = factor(Movement, levels=c("AB Move", "Climate Move"), labels=c("Base", "Climate"))) %>%
    group_by(Recruitment, Movement, hcr, region) %>%
    median_qi(spbio) %>%
    select(Recruitment, Movement, hcr, region, spbio) %>%

    ggplot()+
        geom_col(aes(y=spbio, x=Movement, fill=region))+
        scale_color_manual(values=hcr_colors)+
        facet_grid(Recruitment~hcr)+
        labs(title="Spawning Biomass by Movement Scenario", y="SSB (1000s mt)")+
        custom_theme



#############################################
#############################################
#############################################
#### Landed Catch ###########################
#############################################
#############################################
#############################################
catch_data <- get_landed_catch(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- catch_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )

write_csv(catch_data, file.path("data", "zahneretal_2026_spatialmse_catch.csv"))

catch_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_catch.csv")) %>% mutate(
    Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
    Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
    region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
)


for(r in recruitment_scenarios){
    fname <- paste0("catch_fullutil_", sub(" ", "_", tolower(r)))
    make_catch_plot(
        catch_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

for(r in recruitment_scenarios){
    fname <- paste0("catch_baseutil_", sub(" ", "_", tolower(r)))
    make_catch_plot(
        catch_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)# %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

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


ggsave(filename=file.path(figures_dir, paste0("catch", filetype)), width=width_large, height=height_large, units=c("in"))

catch_props <- catch_data %>% 
    filter(fleet == "Fixed") %>%
    select(-fleet, -catch) %>%
    pivot_wider(names_from=region, values_from=total_catch) %>%
    rowwise() %>%
    mutate(alaska = sum(across(BS:EGOA))) %>%
    mutate(across(BS:EGOA, ~ .x/alaska)) %>%
    filter_times(time_horizon) %>%
    pivot_longer(BS:EGOA, names_to="region", values_to="prop") %>%
    group_by(om, hcr, region) %>%
    median_qi(prop)

catch_props %>% select(om, hcr, region, prop) %>%
    mutate(
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA")),
        om = factor(om, levels=rev(publication_oms))
    ) %>%

    ggplot()+
        geom_bar(aes(y=om, x=prop, fill=region), stat="identity", position="fill")+
        scale_color_manual(values=hcr_colors)+
        facet_wrap(~hcr)+
        custom_theme

catch_data %>% filter_times(time_horizon) %>%
    mutate(Movement = factor(Movement, levels=c("AB Move", "Climate Move"), labels=c("Base", "Climate"))) %>%
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
econ_data <- get_dynamic_economic_value(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
econ_data <- econ_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )

write_csv(econ_data, file.path("data", "zahneretal_2026_spatialmse_econvalue.csv"))

econ_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_econvalue.csv")) %>% 
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    )

for(r in recruitment_scenarios){
    fname <- paste0("economic_value_fullutil_", sub(" ", "_", tolower(r)))
    make_econ_plot(
        econ_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

for(r in recruitment_scenarios){
    fname <- paste0("economic_value_baseutil_", sub(" ", "_", tolower(r)))
    make_econ_plot(
        econ_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)# %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

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

econ_props <- econ_data %>% 
    pivot_wider(names_from=region, values_from=total_value) %>%
    rowwise() %>%
    mutate(alaska = sum(across(BS:EGOA))) %>%
    mutate(across(BS:EGOA, ~ .x/alaska)) %>%
    filter_times(time_horizon) %>%
    pivot_longer(BS:EGOA, names_to="region", values_to="prop") %>%
    group_by(om, hcr, region) %>%
    median_qi(prop)

econ_props %>% select(om, hcr, region, prop) %>%
    mutate(
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA")),
        om = factor(om, levels=rev(publication_oms))
    ) %>%

    ggplot()+
        geom_bar(aes(y=om, x=prop, fill=region), stat="identity", position="fill")+
        scale_color_manual(values=hcr_colors)+
        facet_wrap(~hcr)+
        custom_theme

econ_data %>% filter_times(time_horizon) %>%
    mutate(Movement = factor(Movement, levels=c("AB Move", "Climate Move"), labels=c("Base", "Climate"))) %>%
    group_by(Recruitment, Movement, hcr, region) %>%
    median_qi(total_value) %>%
    select(Recruitment, Movement, hcr, region, total_value) %>%

    ggplot()+
        geom_col(aes(y=total_value, x=Movement, fill=region))+
        scale_color_manual(values=hcr_colors)+
        facet_grid(Recruitment~hcr)+
        labs(title="Total Economic Value by Movement Scenario", y="Total Value")+
        custom_theme


# Plot Recruitment
recruit_data <- get_recruits(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)

recruit_data <- recruit_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    ) %>%
    mutate(region = as.factor(ifelse(is.na(region), "Alaska", as.character(region))))


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

#############################################
#############################################
#############################################
#### ABC, TAC, Expected Landings ############
#############################################
#############################################
#############################################
abc_tac_land <- get_management_quantities(model_runs=NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
abc_tac_land <- abc_tac_land %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )


abc_tac_land <- read_csv(file.path("data", "zahneretal_2026_spatialmse_abctac.csv")) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    )

abc_tac_land <- abc_tac_land %>% filter(fleet == "Trawl", L1 != "tac")

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



#############################################
#############################################
#############################################
#### Average Population Age #################
#############################################
#############################################
#############################################
average_age_data <- get_average_age(model_runs = NULL, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)
average_age_data <- average_age_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )


average_age_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_avgage.csv")) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    )


for(r in recruitment_scenarios){
    fname <- paste0("average_age_fullutil_", sub(" ", "_", tolower(r)))
    make_avgage_plot(
        average_age_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

for(r in recruitment_scenarios){
    fname <- paste0("average_age_baseutil_", sub(" ", "_", tolower(r)))
    make_avgage_plot(
        average_age_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

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

ggsave(file.path(here::here(), "figures", "average_age.jpeg"), width=width_large, height=height_large, units="in")


#############################################
#############################################
#############################################
#### Numbers at Price Grade #################
#############################################
#############################################
#############################################
napg_data <- get_numbers_at_age(model_runs=NULL, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms, pop_or_catch = "catch")
napg_data <- napg_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )

write_csv(napg_data, file.path("data", "zahneretal_2026_spatialmse_napg_cat.csv"))

napg_data <- read_csv(file.path("data", "zahneretal_2026_spatialmse_napg.csv")) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    )

for(r in recruitment_scenarios){
    fname <- paste0("napg_fullutil_", sub(" ", "_", tolower(r)))
    make_napg_plot(
        napg_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

for(r in recruitment_scenarios){
    fname <- paste0("napg_baseutil_", sub(" ", "_", tolower(r)))
    make_napg_plot(
        napg_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

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

ggsave(file.path(here::here(), "figures", "napg.jpeg"), width=width_large, height=height_large, units="in")


###### Performance Metrics
time_horizon <- c(55, 105)
tic()
performance_metrics <- performance_metric_summary(
    model_runs=NULL, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr", "region"),
    hcr_filter=base_util_hcrs,
    om_filter=publication_oms,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "avg_catch_lg", "dynamic_value")
)
toc()
# perf_data <- performance_metrics$perf_data %>% filter(hcr %in% publication_hcrs)
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion Large Catch", "Dynamic Annual Value")
publication_metrics2 = c("Annual\nCatch", "Catch\nAAV", "SSB", "Average\nAge", "Average Catch\nLarge", "Dynamic\nValue")
perf_data2 <- performance_metrics$perf_data %>% 
    # filter(hcr %in% hcrs, name %in% oms) %>%
    filter_hcr_om(publication_hcrs, publication_oms) %>%
    mutate(
            hcr = factor(hcr, levels=publication_hcrs),
            om = factor(om, labels=publication_oms),
            name = factor(name, levels=publication_metrics, labels=publication_metrics2)
    ) %>%
    filter(.width == 0.5) %>%
    group_by(region, om, name) %>%
    scale_and_rank("median")

perf_data2 %>% select(om, hcr, region, name, scaled) %>% arrange(om, region, name) %>%
    ggplot(aes(x=name, y=scaled, color=hcr, fill=hcr, group=hcr))+
        geom_point(size=3)+
        geom_line()+
        geom_polygon(alpha=0)+
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))+
        scale_color_manual(values=hcr_colors)+
        ggiraphExtra::coord_radar()+
        facet_grid(region~om)+
        custom_theme+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
ggsave(file.path(here::here(), "figures", "performance_metrics_radar.png"), width=25, height=25, units="in")


##### Resiliency Metrics
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

