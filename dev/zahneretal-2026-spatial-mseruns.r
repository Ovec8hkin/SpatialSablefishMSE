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
nsims <- 100
seed_list <- sample(1:(10*nsims), nsims)  # Draw 10 random seeds

mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 54
mse_options$recruitment_start_year <- 54
mse_options$n_proj_years <- 50
mse_options$run_estimation <- TRUE

mse_options_list <- listN(mse_options)

om_list <- listN(om_agemove_bhrec, om_agemove_regimerec, om_agemove_crashrec)#, om_agemove_crashrec)#, om_agemove_regimerec, om_climatemove_regimerec, om_agemove_crashrec, om_climatemove_crashrec)
hcr_list <- listN(mp_f40, mp_f40hybrid)#, mp_dynamicB0)#, mp_f40_fullutil, mp_f50, mp_f50_fullutil, mp_b30f50, mp_b30f50_fullutil, mp_f40hybrid, mp_f40hybrid_fullutil)

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
extra_columns <- expand.grid(hcr=lapply(hcr_list, function(x) x$name), om=lapply(om_list, function(x) x$name))


# Data Processing
filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures")
width_small <- 12
height_small <- 8

width_large <- 16
height_large <- 16

width_medium <- 12
height_medium <- 12

base_util_hcrs <- c("F40", "F50", "F50/B30", "F40 Hybrid")
full_util_hcrs <- paste0(base_util_hcrs, " | Full Utilization")

publication_hcrs <- c("F40", "F40 | Full Utilization", "F50", "F50 | Full Utilization", "F50/B30", "F50/B30 | Full Utilization", "F40 Hybrid", "F40 Hybrid | Full Utilization")
publication_oms <- c("BH Recruit | AB Move", "BH Recruit | Climate Move", "Regime Recruit | AB Move", "Regime Recruit | Climate Move", "Crash Recruit | AB Move", "Crash Recruit | Climate Move")
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years with Low SSB")

recruitment_scenarios <- c("BH Recruit", "Regime Recruit", "Crash Recruit")

extra_columns = expand.grid(
    om=c("BH Recruit | AB Move", "Regime Recruit | AB Move", "Crash Recruit | AB Move"),
    hcr=c("F40Test", "F40 Hybrid")
)

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 105)

hcr_colors <- c("black", "#991c1c", "#001180", "#256c15", "#BB6A00", "#7500B0", "#5CABA3", "#BD4E98")
#names(hcr_colors) <- publication_hcrs

hcr_colors = c("black", "#256c15", "#0763d3", "#BD4E98")

hcr_filter <- c("F40", "F40 Hybrid")
om_filter <- c("BH Recruit | AB Move", "Regime Recruit | AB Move", "Crash Recruit | AB Move")

###### Plot Diagnostics
# plot_em_diagnostics(
#     model_runs, 
#     extra_columns,
#     publication_hcrs,
#     publication_oms,
#     spinup_years = 54,
#     n_proj_years = 50,
#     simulation_year = 4,
#     simulation_number = 1
# )

convergence <- check_em_convergence_diagnostics(
    model_runs=NULL,
    extra_columns,
    publication_hcrs,
    publication_oms
)

unconverged <- convergence %>% filter(!is.na(failed))

# convergence_summ <- convergence %>% 
#     mutate(
#         sim=rep(seed_list, each=51),
#         year = rep(54:(54+50), 15)
#     )

# convergence_summ %>% group_by(year) %>% summarise(c = sum(converged)/n()) %>% print(n=100)

plot_diags(
    hcr = "F40Test",
    om = "Regime Recruit | AB Move",
    simulation_seed = 1602,
    simulation_year = 1,
    om_list, 
    seed_list
)

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

ssb_data %>% select(time, sim, L1, om, Recruitment, Movement, hcr, region, spbio) %>%
    filter_times(time_horizon) %>%
    pivot_wider(names_from=L1, values_from=spbio) %>%
    group_by(time, sim, om, Recruitment, Movement, hcr) %>%
    mutate(naa = sum(naa, na.rm=TRUE)) %>%
    filter(region == "Alaska") %>%
    ungroup() %>%
    mutate(
        bias = (naa_est - naa)/naa
    ) %>%
    # group_by(time, om, Recruitment, Movement, hcr) %>%
    # median_qi(bias, .width=c(0.50, 0.80)) %>%

    ggplot()+
        geom_boxplot(aes(x=om, y=bias))+
        geom_hline(yintercept=0, linetype="dashed")+
        facet_wrap(~hcr)+
        scale_y_continuous(breaks=seq(-0.1, 0.2, 0.02))+
        custom_theme

    # ggplot()+
    #     geom_lineribbon(aes(x=time, y=bias, ymin=.lower, ymax=.upper, group=hcr))+
    #     geom_hline(yintercept=0, linetype="dashed")+
    #     scale_fill_brewer(palette="Blues")+
    #     facet_wrap(~om)+
    #     custom_theme



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
            scale_y_continuous(limits=c(0, 75), breaks=seq(0, 75, 25)),
            scale_y_continuous(limits=c(0, 75), breaks=seq(0, 75, 25)),
            scale_y_continuous(limits=c(0, 40), breaks=seq(0, 40, 15)),
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
            v1="hcr", v2="om", v3=NA, show_est = TRUE,
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

##### SSB Timeseries Relative to the F40 Harvest Control Rule
reg_ssb_f40rel_plots <- plot_ssb(
        ssb_data, 
        v1="hcr", v2="region", v3="om", show_est = FALSE,
        common_trajectory=common_trajectory, time_horizon = time_horizon,
        relative="F40"
    )+
    labs(title="Regional Spawning Biomass Relative to F40", y="Relative SSB")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    scale_y_continuous(limits=c(0.5, 2))+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

agg_ssb_f40rel_plots <- plot_ssb(
        ssb_data, 
        v1="hcr", v2="om", v3=NA, show_est = TRUE,
        common_trajectory=common_trajectory, time_horizon = time_horizon,
        relative="F40"
    )+
    scale_y_continuous("Relative SSB", limits=c(0.5, 2))+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_ssb_f40rel_plots / agg_ssb_f40rel_plots +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")

ggsave(filename=file.path(figures_dir, paste0("ssb_f40rel", filetype)), width=width_large, height=height_large, units="in")

##### SSB Timeseries as difference between climate and age-based movement
#########################################################################

ab_ssb <- ssb_data %>% filter(Movement == "AB Move")
cl_ssb <- ssb_data %>% filter(Movement == "Climate Move")

ssb_diff <- ab_ssb %>% 
    left_join(cl_ssb, by=c("time", "sim", "Recruitment", "hcr", "region")) %>%
    mutate(
        spbio = spbio.y - spbio.x,
        biomass = biomass.y - biomass.x
    ) %>%
    select(time, sim, L1=L1.x, Recruitment, hcr, region, spbio, biomass)


depletion_diff_plots <- plot_ssb(ssb_diff, v1="hcr", v2="region", v3="Recruitment", common_trajectory=common_trajectory, show_est = FALSE, scales="free")+
    labs(title="Regional Depletion (SSB/SSB0)", y="Depletion")+
    facet_grid(region~Recruitment, scales="free_y")+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 25)),
            scale_y_continuous(limits=c(0, 75), breaks=seq(0, 75, 25)),
            scale_y_continuous(limits=c(0, 25), breaks=seq(0, 25, 5)),
            scale_y_continuous(limits=c(-75, 0), breaks=seq(-75, 0, 25)),
            scale_y_continuous(limits=c(-75, 0), breaks=seq(-75, 0, 25))
        )
    )+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

ssb_diff_agg_plots <- plot_ssb(ssb_diff, v1="hcr", v2="Recruitment", show_est=FALSE, common_trajectory=common_trajectory, scales="free_y")+
    scale_y_continuous(limits=c(-5, 25))+
    facet_grid(region~Recruitment, scales="free_y")+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

# ssb_agg_plots + depletion_plots + 
depletion_diff_plots / ssb_diff_agg_plots +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")


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

### Landings Timeseries Relative to the F40 Harvest Control Rule
reg_catch_f40rel_plots <- plot_landed_catch(catch_data, v1="hcr", v2="om", v3="region", by_fleet=FALSE, common_trajectory = common_trajectory, time_horizon=time_horizon, relative="F40")+
    labs(title="Regional Landings Relative to F40", y="Relative Catch")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    scale_y_continuous(limits=c(0, 2.5))+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

agg_catch_f40rel_plot <- plot_landed_catch(
        catch_data, 
        v1="hcr", v2="om", by_fleet=FALSE, 
        common_trajectory = common_trajectory, time_horizon=time_horizon,
        relative="F40"
    )+
    scale_y_continuous(limits=c(0.25, 1.5), name="Relative Catch")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_catch_f40rel_plots / agg_catch_f40rel_plot +  
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")

ggsave(filename=file.path(figures_dir, paste0("catch_f40rel", filetype)), width=width_large, height=height_large, units="in")


ab_catch <- catch_data %>% filter(Movement == "AB Move")
cl_catch <- catch_data %>% filter(Movement == "Climate Move")

catch_diff <- ab_catch %>% 
    left_join(cl_catch, by=c("time", "sim", "fleet", "Recruitment", "hcr", "region")) %>%
    mutate(
        total_catch = total_catch.y - total_catch.x,
        catch = catch.y - catch.x
    ) %>%
    select(time, sim, L1=L1.x, fleet, Recruitment, hcr, region, catch, total_catch)


reg_catch_diff_plots <- plot_landed_catch(catch_diff, v1="hcr", v2="region", v3="Recruitment", common_trajectory=common_trajectory)+
    labs(title="Regional Landings", y="Catch (mt)")+
    facet_grid(region~Recruitment, scales="free_y")+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 5), breaks=seq(0, 5, 1)),
            scale_y_continuous(limits=c(0, 5), breaks=seq(0, 5, 1)),
            scale_y_continuous(limits=c(0, 5), breaks=seq(0, 5, 1)),
            scale_y_continuous(limits=c(-5, 0), breaks=seq(-5, 0, 1)),
            scale_y_continuous(limits=c(-5, 0), breaks=seq(-5, 0, 1))
        )
    )+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

catch_diff_agg_plots <- plot_landed_catch(catch_diff, v1="hcr", v2="Recruitment", common_trajectory=common_trajectory)+
    scale_y_continuous(limits=c(-5, 2))+
    facet_grid(region~Recruitment, scales="free_y")+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_catch_diff_plots / catch_diff_agg_plots +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")





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
econ_data <- get_dynamic_economic_value(model_runs=NULL, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)
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
                scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
                scale_y_continuous(limits=c(0, 3), breaks=seq(0, 2, 0.5)),
                scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
                scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 2)),
                scale_y_continuous(limits=c(0, 6), breaks=seq(0, 6, 2))
            )
        )+
        theme(
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )

    econ_agg_plot <- plot_econ_value(econ_data, v1="hcr", v2="om", common_trajectory = common_trajectory)+
        scale_y_continuous(limits=c(0, 15))+
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

##### Economic Value Timeseries Relative to F40 Harvest Control Rule

reg_econ_f40rel_plots <- plot_econ_value(
        econ_data, 
        v1="hcr", v2="om", v3="region", 
        common_trajectory = common_trajectory, time_horizon=time_horizon,
        relative="F40"
    )+
    labs(title="Regional Economic Value Relative to F40", y="Relative Economic Value")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    scale_y_continuous(limits=c(0.5, 2))+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

agg_econ_f40rel_plot <- plot_econ_value(
        econ_data, 
        v1="hcr", v2="om", 
        common_trajectory = common_trajectory, time_horizon = time_horizon,
        relative="F40"
    )+
    scale_y_continuous(limits=c(0.5, 1.25))+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_econ_f40rel_plots / agg_econ_f40rel_plot +  
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5))+
    plot_annotation(caption="Economic value currently calculated for Fixed gear fleet only") & 
    theme(legend.position = "bottom")

ggsave(filename=file.path(figures_dir, paste0("economic_value_f40rel", filetype)), width=width_large, height=height_large, units="in")




ab_econ <- econ_data %>% filter(Movement == "AB Move")
cl_econ <- econ_data %>% filter(Movement == "Climate Move")

econ_diff <- ab_econ %>% 
    left_join(cl_econ, by=c("time", "sim", "Recruitment", "hcr", "region")) %>%
    mutate(
        total_value = total_value.y - total_value.x,
    ) %>%
    select(time, sim, Recruitment, hcr, region, total_value)


reg_econ_diff_plots <- plot_econ_value(econ_diff, v1="hcr", v2="region", v3="Recruitment", common_trajectory=common_trajectory)+
    labs(title="Regional Economic Value", y="Value")+
    facet_grid(region~Recruitment, scales="free_y")+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
            scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
            scale_y_continuous(limits=c(0, 2), breaks=seq(0, 2, 0.5)),
            scale_y_continuous(limits=c(-3, 0), breaks=seq(-3, 0, 1)),
            scale_y_continuous(limits=c(-3, 0), breaks=seq(-3, 0, 1))
        )
    )+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

econ_diff_agg_plots <- plot_econ_value(econ_diff, v1="hcr", v2="Recruitment", common_trajectory=common_trajectory)+
    scale_y_continuous(limits=c(-3, 1))+
    facet_grid(region~Recruitment, scales="free_y")+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_econ_diff_plots / econ_diff_agg_plots +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")

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
recruit_data <- get_recruits(model_runs, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)

reg_rec_plots <- plot_recruitment(recruit_data, v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)

recruit_agg_plot <- plot_recruitment(recruit_data, v1="hcr", v2="om", common_trajectory = common_trajectory)+
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



#############################################
#############################################
#############################################
#### ABC, TAC, Expected Landings ############
#############################################
#############################################
#############################################
abc_tac_land <- get_management_quantities(model_runs=NULL, extra_columns, hcr_filter=c("F40Test", "F40 Hybrid"), om_filter=c("BH Recruit | AB Move", "Regime Recruit | AB Move", "Crash Recruit | AB Move"))
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
    # scale_color_manual(values=c("red", "blue", "green", "purple"))+
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

###### Average Age Tiemseris Relative to F40 Harvest Control Rule

reg_avgage_f40rel_plots <- plot_average_age(
        average_age_data, 
        v1="hcr", v2="om", v3="region", 
        common_trajectory = common_trajectory, time_horizon=time_horizon,
        relative="F40"
    )+
    labs(title="Average Population Age Relative to F40", y="Relative Age")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    scale_y_continuous(limits=c(0.9, 1.2))+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

agg_avgage_f40rel_plot <- plot_average_age(
        average_age_data, 
        v1="hcr", v2="om", v3=NA, 
        common_trajectory = common_trajectory, time_horizon = time_horizon,
        relative="F40"
    )+
    # scale_color_manual(values=c("red", "blue", "green", "purple"))+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    guides(color="none", shape="none")+
    labs(y="Relative Age")+
    scale_y_continuous(limits=c(0.9, 1.2), breaks=seq(0.9, 1.2, 0.1))+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_avgage_f40rel_plots / agg_avgage_f40rel_plot +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")
ggsave(file.path(here::here(), "figures", "average_age_f40rel.jpeg"), width=width_large, height=height_large, units="in")



#############################################
#############################################
#############################################
#### Numbers at Price Grade #################
#############################################
#############################################
#############################################
napg_data <- get_numbers_at_age(model_runs=NULL, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)
napg_data <- napg_data %>% separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
    )

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
                scale_y_continuous(limits=c(0, 15)),
                scale_y_continuous(limits=c(0, 15)),
                scale_y_continuous(limits=c(0, 10)),
                scale_y_continuous(limits=c(0, 25)),
                scale_y_continuous(limits=c(0, 25))
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
        scale_y_continuous(limits=c(0, 45))+
        theme(
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )

    reg_highpg_plots / agg_highpg_plot +
        plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
        theme(legend.position = "bottom")
}

ggsave(file.path(here::here(), "figures", "napg.jpeg"), width=width_large, height=height_large, units="in")

### Numbers at High Price Grade Timeseries Relative to F40 Harvest Control Rule

reg_highpg_f40rel_plots <- plot_numbers_at_pricegrade(
        napg_data, 
        v1="hcr", v2="om", v3="region", price_grade="Grade 7+ (15+yo)", 
        common_trajectory = common_trajectory, time_horizon=time_horizon,
        relative="F40"
    )+
    labs(title="Average Number Grade 7+ Individuals Relative to F40", y="Relative Individuals")+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    scale_y_continuous(limits=c(0.5, 2.5))+
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
    )

agg_highpg_f40rel_plot <- plot_numbers_at_pricegrade(
        napg_data, 
        v1="hcr", v2="om", v3=NA, price_grade="Grade 7+ (15+yo)", 
        common_trajectory = common_trajectory, time_horizon=time_horizon,
        relative="F40"
    )+
    ggh4x::facet_nested(
        rows=vars(region), 
        cols=vars(Recruitment, Movement), 
        scales="free_y"
    )+
    guides(color="none", shape="none")+
    labs(y="Relative Individuals")+
    scale_y_continuous(limits=c(0.5, 2.25))+
    theme(
        strip.background.x = element_blank(),
        strip.text.x = element_blank()
    )

reg_highpg_f40rel_plots / agg_highpg_f40rel_plot +
    plot_layout(ncol=1, guides="collect", heights=c(1, 0.5)) & 
    theme(legend.position = "bottom")

ggsave(file.path(here::here(), "figures", "napg_f40rel.jpeg"), width=width_large, height=height_large, units="in")


















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


full_performance_metrics <- performance_metric_summary(
    model_runs_saved, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=publication_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr", "region"),
    summary_out = FALSE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "avg_catch_lg", "dynamic_value") 
)

om_aggregated_performance <- full_performance_metrics$perf_data %>%
    separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
    pivot_longer(7:12, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "catch", "dyn_annual_value"), 
            labels=publication_metrics
        )
    ) %>%
    pivot_wider(names_from="name", values_from="value")
    
recruitment_aggregated_performance <- om_aggregated_performance %>%   
    group_by(sim, hcr, Recruitment, region) %>%
    summarise(
        across(3:8, \(x) mean(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(5:10, names_to="name", values_to="value")

movement_aggregated_performance <- om_aggregated_performance %>%   
    group_by(sim, hcr, Movement, region) %>%
    summarise(
        across(3:8, \(x) mean(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(5:10, names_to="name", values_to="value")



movement_aggregated_performance %>% 
    group_by(Movement, hcr, region, name) %>%
    median_qi(value, .width=c(0.50)) %>%
    group_by(Movement, region, name) %>%
    mutate(
        scaled = case_when(
            name == "Catch\nAAV" ~ (inf_max(value)-value)/(inf_max(value)-min(value)),
            TRUE ~ value/inf_max(value)
        )
    ) %>%
    select(Movement, hcr, region, name, scaled) %>%
    # pivot_wider(names_from=name, values_from=scaled) %>%

    ggplot(aes(x=name, y=scaled, color=Movement, fill=Movement, group=Movement))+
        geom_point(size=3)+
        geom_line()+
        geom_polygon(alpha=0)+
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))+
        ggiraphExtra::coord_radar()+
        facet_grid(hcr ~ region)+
        custom_theme+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )

topsis_movement <- compute_topsis(
    # movement_aggregated_performance %>% group_by(Movement, hcr, region, name) %>% median_qi(value, .width=c(0.50)) %>% select(Movement, hcr, region, name, value),
    # perf_data2 %>% filter(.width==0.50) %>% select(om, hcr, region, name, value=median),
    movement_aggregated_performance,# %>% ungroup() %>% pivot_longer(7:ncol(.), names_to="name", values_to="value") %>% select(sim, Movement, hcr, region, name, value), 
    topsis_splits = c("sim", "region", "Movement"),
    topsis_weights = c(0.25, 0.125, 0.25, 0.1, 0.125, 0.15),
    topsis_minmax = c("max", "min", "max", "max", "max", "max")
) %>% arrange(region) %>%
    pivot_longer(4:8, names_to="hcr", values_to="value") %>%
    group_by(region, Movement, hcr) %>%
    median_qi(value, .width=interval_widths) %>%
print(n=100)


ggplot(topsis_movement) + 
    geom_pointinterval(aes(x=value, xmin=.lower, xmax=.upper, y=hcr, color=hcr))+
    facet_grid(Movement~region)+
    scale_x_continuous(name="TOPSIS Score", limits=c(0, 1))+
    custom_theme


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
    om_filter=c("Immediate Crash Recruitment")
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
    om_filter=c("Immediate Crash Recruitment")
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

####################################
#### THESE DONT WORK YET!!!!!!! ####
####################################
# Reference Points
f40_tseries <- get_rp_timeseries(model_runs, extra_columns, hcr_filter=hcrs, om_filter=oms, ref_pt="Fref", spr_target=0.40)


b40s <- get_b40_timeseries(model_runs, extra_columns, publication_hcrs, publication_oms)
ggplot(b40s)+
    geom_lineribbon(aes(x=time, y=B40, ymin=.lower, ymax=.upper))+
    scale_fill_brewer(palette="Blues")+
    scale_y_continuous(limits=c(0, 200))+
    facet_wrap(~om)+
    custom_theme

# Plot recruitment
rec_data <- get_recruits(model_runs, extra_columns2, hcr_filter=hcr_names, om_filter=om_names)
plot_recruitment(rec_data, v1="hcr", v2="om", show_est = FALSE)+custom_theme



hcr_fs <- bind_mse_outputs(model_runs, "hcr_f", extra_columns)

hcr_fs_tseries <- hcr_fs %>% filter_times(time_horizon) %>% filter_hcr_om(publication_hcrs, publication_oms) %>%
    as_tibble() %>%
    select(-Var2, -Var3, -region, -L1) %>%
    group_by(time, om, hcr) %>%
    median_qi(value, .width=interval_widths) 

ggplot(hcr_fs_tseries, aes(x=time, y=value, color=hcr))+
    geom_line()+
    facet_wrap(~om)+
    custom_theme
