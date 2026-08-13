rm(list=ls())
    
library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SPoRC)
library(tictoc)
library(doParallel)
library(patchwork)
library(afscOM) # may work but not certain
library(devtools)

sablefishMSE_dir <- here::here()
lapply(list.files("R", full.names = TRUE), source)

source("dev/hcrs.R")

nyears <- 110
sable_om <- readRDS("data/spatial_sablefish_om.RDS") # Read this saved OM from a file
assessment <- dget(file.path(here::here(), "data/spatial_sablefsh_inputs.rdat"))
hist_recruits <- t(assessment$recruitment)

ricker_recruit <- function(a, b, sigR, seed){
    set.seed(seed)
    function(ssb, y){
        r <- a * ssb * exp(-b * ssb)
        return(r)
    }
}

###########################
###########################
####### Create new OM using Ricker stock recruitment functions
###########################
###########################
om_ricker_recruit <- sable_om
om_ricker_recruit$name <- "Ricker Recruitment" # Name of Operating Model Object
om_ricker_recruit$recruitment$func <- ricker_recruit # Recruitment function to generate new recruits from
# Required parameters for the Ricker recruitment function
om_ricker_recruit$recruitment$pars <- list(
    a=1, 
    b=0.001, 
    sigR= 0.7
)
om_ricker_recruit$recruitment$apportionment$func <- resample_recruit_apportionment # Apportionment recruits by random resampling
om_ricker_recruit$recruitment$apportionment$pars <- list(
    hist_recruits = hist_recruits
)

# For effect, turn off generating RPWs for the GOA trawl survey (Fleet 4)
om_ricker_recruit$model_options$obs_pars$is_survey <- c(0, 0, 1, 1)
om_ricker_recruit$model_options$obs_pars$rpw <- c(0, 0, 1, 0)
om_ricker_recruit$model_options$obs_pars$rpw_cv <- c(0, 0, 0.1, 0)
om_ricker_recruit$model_options$obs_pars$acs <- c(1, 1, 1, 1)
om_ricker_recruit$model_options$obs_pars$ac_samps <- c(200, 200, 200, 200)


###########################
###########################
####### Create new OM with regionally varying WAA
###########################
###########################
om_regional_dems <- sable_om
om_regional_dems$name <- "Regional Demographics"
om_regional_dems$recruitment <- om_ricker_recruit$recruitment

custom_waa <- array(
    c(
        seq(1, 4, length.out=30),   # BS Females
        seq(1, 2.5, length.out=30), # BS Males
        seq(1, 4, length.out=30),   # AI Females
        seq(1, 2.5, length.out=30), # AI Males
        seq(1, 5, length.out=30),   # WGOA Females
        seq(1, 3.5, length.out=30), # WGOA Males
        seq(1, 6, length.out=30),   # CGOA Females
        seq(1, 3.5, length.out=30), # CGOA Males
        seq(1, 8, length.out=30),   # EGOA Females
        seq(1, 5.5, length.out=30) # EGOA Males
    ),
    dim=c(30, 2,5),
    dimnames=list(
        "age"=2:31,
        "sex"=c("F", "M"),
        "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA")
    )
)

waa <- afscOM::generate_param_matrix(
    custom_waa, 
    dimension_names = list(
        "time"=55:200,
        "age"=2:31,
        "sex"=c("F", "M"),
        "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA")
    ),
    by=c("age", "sex", "region")
)

om_regional_dems$dem_params$waa[55:200,,,] <- waa



###########################
###########################
####### Create new MP Object
###########################
###########################
mp_base <- setup_mp_options()

# Setup a new HCR using F reference level of SPR50,
# biomass reference level of SPR30, with a harvest
# cap of 15,000 mt, and a bidirectional 15% stability
# constraint.
mp_newhcr <- mp_base
mp_newhcr$name <- "New HCR"

mp_newhcr$hcr <- list(
    func = tier3, # Standard rectilinear HCR function, see hcrs.R
    extra_pars = NA, # No additional components requires for parameterizing HCR function
    extra_options = list(
        max_stability = array(c(0.15, 0.15), dim=c(1, 2)),
        harvest_cap = 15 # In thousands of metric tons
    ),
    units = "F"
)
mp_newhcr$ref_points$spr_target <- c(0.50, 0.30) # F ref point, then B rep point

mp_newhcr$management$abc_tac_regflt_reduction <- matrix(
    data = 0.95,
    nrow=5,
    ncol=2,
    dimnames=list("region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl"))
)




om_list <- listN(om_ricker_recruit, om_regional_dems)
hcr_list <- listN(mp_f40, mp_newhcr)

mse_options <- setup_mse_options()
mse_options$n_proj_years <- 15
mse_options$run_estimation <- FALSE
mse_options$mse_verbose <- TRUE

mse_options_list <- listN(mse_options)

set.seed(1120)
nsims <- 10
seed_list <- sample(1:(10*nsims), nsims)  # Draw 10 random seeds

model_runs <- run_mse_multiple(
    om_list,
    hcr_list,
    seed_list,
    mse_options_list,
    diagnostics=FALSE,
    out.dir=file.path(here::here(), "data", "test"),
    save=TRUE,
    overwrite = TRUE
)

hcr_filter <- c("F40", "New HCR")
om_filter <- c(unlist(lapply(om_list, \(x) x$name)))

extra_columns <- expand.grid(
    om = om_filter,
    hcr = hcr_filter
)

interval_widths <- c(0.50, 0.95)
common_trajectory <- 54
time_horizon <- c(55, 70)

ssb <- get_ssb_biomass(
    model_runs=NULL, 
    extra_columns=extra_columns,
    dem_params=sable_om$dem_params,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    dir=file.path("data", "test")
)

hcr_colors <- c("black", "blue")
plot_ssb(ssb %>% filter(time > 45), v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)+
    ggh4x::facet_nested(
        rows=vars(om), 
        cols=vars(region), 
        scales="free"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 650)),
            scale_y_continuous(limits=c(0, 300))
        )
    )

catch <- get_landed_catch(
    model_runs=NULL, 
    extra_columns=extra_columns,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    dir=file.path("data", "test")
)
plot_landed_catch(catch %>% filter(time > 45), v1="hcr", v2="om", v3="region", common_trajectory = common_trajectory)+
    ggh4x::facet_nested(
        rows=vars(om), 
        cols=vars(region), 
        scales="free"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 70)),
            scale_y_continuous(limits=c(0, 70))
        )
    )

perf_metrics <- performance_metric_summary(
    model_runs=NULL, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = NULL,
    interval_widths=c(0.50, 0.80),
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr", "region"),
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb"),
    dir=file.path("data", "test")
)

perf_data <- perf_metrics$perf_data %>%
    # filter(.width == 0.5, region != "Alaska") %>%
    group_by(om, region, name) %>%
    scale_and_rank("median") %>%
    mutate(
            hcr = factor(hcr, levels=hcr_filter),
            om = factor(om, levels=om_filter)
            # name = factor(name, levels=c("Annual Catch", "Catch AAV", "SSB"), labels=publication_metric_names)
    )

perf_data %>% 
    select(om, hcr, region, name, scaled) %>% 
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
            cols=vars(om), 
            rows=vars(region)
        )+
        custom_theme+
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            # panel.spacing.x = grid::unit(rep_len(c(0.1, 0.5), 5), "in")
            # panel.spacing.y = unit(c(0.1, 0.1, 0.1, 0.1, 0.5), "in")        
        )
