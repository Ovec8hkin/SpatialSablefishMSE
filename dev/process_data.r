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


#### Extra Plotting Functions
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
            v1="hcr", v2="om", v3=NA, show_est = FALSE,
            common_trajectory=common_trajectory
        )+
        scale_y_continuous(limits=c(0, 300))+
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
    om=publication_oms,
    hcr=publication_hcrs
)

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 105)

#               F40        F50          F50/B30      F40 Hybrid
hcr_colors = c("black", "#256c15", "#0763d3", "#BD4E98")

# Data files
ssb_data_file       <- file.path("data", "zahneretal_2026_spatialmse_SSB_biomass.csv")
catch_data_file     <- file.path("data", "zahneretal_2026_spatialmse_catch.csv")
econ_data_file      <- file.path("data", "zahneretal_2026_spatialmse_econvalue.csv")
avgage_data_file    <- file.path("data", "zahneretal_2026_spatialmse_avgage.csv")
napg_data_file      <- file.path("data", "zahneretal_2026_spatialmse_napg.csv")
abc_data_file       <- file.path("data", "zahneretal_2026_spatialmse_abctac.csv")

#############################################
#############################################
#############################################
#### Spawning Biomass + Depletion ###########
#############################################
#############################################
#############################################
ssb_data <- read_csv(ssb_data_file) %>% set_factor_levels

# Full Utilization HCRs
for(r in recruitment_scenarios){
    fname <- paste0("ssb_depletion_fullutil_", sub(" ", "_", tolower(r)))
    make_ssb_plot(
        ssb_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_ssb_plot(
    ssb_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("ssb_depletion_fullutil", filetype)), width=width_large, height=height_medium, units=c("in"))

# Base Utilization HCRs
for(r in recruitment_scenarios){
    fname <- paste0("ssb_depletion_baseutil_", sub(" ", "_", tolower(r)))
    make_ssb_plot(
        ssb_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_ssb_plot(
    ssb_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("ssb_depletion_baseutil", filetype)), width=width_large, height=height_medium, units=c("in"))




#############################################
#############################################
#############################################
#### Landed Catch ###########################
#############################################
#############################################
#############################################
catch_data <- read_csv(catch_data_file) %>% set_factor_levels()

for(r in recruitment_scenarios){
    fname <- paste0("catch_fullutil_", sub(" ", "_", tolower(r)))
    make_catch_plot(
        catch_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_catch_plot(
    catch_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("catch_fullutil", filetype)), width=width_large, height=height_medium, units=c("in"))


for(r in recruitment_scenarios){
    fname <- paste0("catch_baseutil_", sub(" ", "_", tolower(r)))
    make_catch_plot(
        catch_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_catch_plot(
    catch_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("catch_baseutil", filetype)), width=width_large, height=height_medium, units=c("in"))




#############################################
#############################################
#############################################
#### Dynamic Economic Value #################
#############################################
#############################################
#############################################
econ_data <- read_csv(econ_data_file) %>% set_factor_levels

for(r in recruitment_scenarios){
    fname <- paste0("economic_value_fullutil_", sub(" ", "_", tolower(r)))
    make_econ_plot(
        econ_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_econ_plot(
    econ_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("catch_fullutil", filetype)), width=width_large, height=height_medium, units=c("in"))


for(r in recruitment_scenarios){
    fname <- paste0("economic_value_baseutil_", sub(" ", "_", tolower(r)))
    make_econ_plot(
        econ_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_econ_plot(
    econ_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("econ_baseutil", filetype)), width=width_large, height=height_medium, units=c("in"))



econ_props <- econ_data %>% rename(value="total_value") %>% 
    filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>%
    filter_times(time_horizon) %>%
    group_by(time, sim, om, hcr) %>%
    mutate(
        prop = value/sum(value)
    ) %>%
    group_by(time, Recruitment, Movement, hcr, region) %>%
    median_qi(prop, .width=c(0.50))

    ggplot(econ_props)+
        geom_bar(aes(x=time, y=prop, fill=region), stat="identity", position="fill")+
        geom_hline(yintercept=c(0.25, 0.50, 0.75), linetype="dashed")+
        ggh4x::facet_nested(
            rows = vars(hcr),
            cols = vars(Recruitment, Movement)
        )+
        scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))+
        scale_fill_manual(values=c("red", "blue", "yellow", "darkgreen", "purple"))+
        labs(x = "Time", y="Proportion of Total Economic Value", fill="Region", title = "Regional Proportion of Total Economic Value")+
        custom_theme

ggsave(file.path(figures_dir, "regprop_plots", paste0("econvalue_baseutil_prop", filetype)), width=width_large, height=height_large, units="in")




#############################################
#############################################
#############################################
#### Numbers at Price Grade #################
#############################################
#############################################
#############################################
napg_data <- read_csv(napg_data_file) %>% set_factor_levels


for(r in recruitment_scenarios){
    fname <- paste0("napg_fullutil_", sub(" ", "_", tolower(r)))
    make_napg_plot(
        napg_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_napg_plot(
    napg_data %>% filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("napg_fullutil", filetype)), width=width_large, height=height_medium, units=c("in"))


for(r in recruitment_scenarios){
    fname <- paste0("napg_baseutil_", sub(" ", "_", tolower(r)))
    make_napg_plot(
        napg_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms) %>% filter(Recruitment == r)
    )
    ggsave(filename=file.path(figures_dir, paste0(fname, filetype)), width=width_medium, height=height_medium, units=c("in"))
}

make_napg_plot(
    napg_data %>% filter_hcr_om(hcrs=base_util_hcrs, oms=publication_oms)
)
ggsave(filename=file.path(figures_dir, paste0("napg_baseutil", filetype)), width=width_large, height=height_medium, units=c("in"))



napg_props <- napg_data %>% filter(grepl("Grade 7+", class)) %>%
    filter_hcr_om(hcrs=full_util_hcrs, oms=publication_oms) %>%
    filter_times(time_horizon) %>%
    group_by(time, sim, om, hcr) %>%
    mutate(
        prop = value/sum(value)
    ) %>%
    group_by(time, Recruitment, Movement, hcr, region) %>%
    median_qi(value, .width=c(0.50))

    ggplot(napg_props)+
        geom_col(aes(x=time, y=value, fill=region))+
        # geom_hline(yintercept=c(0.25, 0.50, 0.75), linetype="dashed")+
        ggh4x::facet_nested(
            rows = vars(hcr),
            cols = vars(Recruitment, Movement)
        )+
        scale_y_continuous(limits=c(0, 45), breaks=seq(0, 40, 10))+
        scale_fill_manual(values=c("red", "blue", "yellow", "darkgreen", "purple"))+
        labs(y="Proportion of Total Grade 7+ Individuals", x="Time", fill="Region", title="Regional Proportion of High Price Grade Individuals")+
        custom_theme

ggsave(file.path(figures_dir, "regprop_plots", paste0("napg_fullutil", filetype)), width=width_large, height=height_large, units="in")



set_factor_levels <- function(data){
    data %>% mutate(
        Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
        Movement = factor(Movement, levels=c("AB Move", "Climate Move")),
        region = factor(region, levels=c("BS", "AI", "WGOA", "CGOA", "EGOA"))
    )
}
