rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(tictoc)
library(abind)

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"

devtools::load_all(afscOM_dir)
# devtools::load_all(sablefishMSE_dir)

lapply(list.files("R", full.names = TRUE), source)

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons

spatial_sablefish_data <- readRDS("~/Desktop/Projects/SabieTMB/output/Spatial Assessment/data.RDS")
spatial_sablefish_estimates <- readRDS("~/Desktop/Projects/SabieTMB/output/Spatial Assessment/mle_rep.RDS")

nyears <- length(spatial_sablefish_data$years)
nages <- length(spatial_sablefish_data$ages)
nsexes <- 2
nregions <- 5
nfleets <- 2
nsurveys <- 2

dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 2:(nages+1),
    "sex"  = c("F", "M"),
    "region" = c("BS", "AI", "WGOA", "CGOA", "EGOA"),
    "fleet" = c("Fixed", "Trawl")
)

model_params <- afscOM::set_model_params(nyears, nages, nsexes, nregions, nfleets)
model_options <- afscOM::setup_model_options(model_params)

# 1. Generate demographic parameter matrices
mort_matrix <- afscOM::generate_param_matrix(
    aperm(spatial_sablefish_estimates$natmort, perm=c(2, 3, 4, 1)),
    by=c("time", "age", "sex", "region"),
    dimension_names = dimension_names
)
waa_matrix <- afscOM::generate_param_matrix(
    aperm(spatial_sablefish_data$WAA, perm=c(2, 3, 4, 1)),
    by=c("time", "age", "sex", "region"),
    dimension_names = dimension_names
)
mat_matrix <- afscOM::generate_param_matrix(
    aperm(spatial_sablefish_data$MatAA, perm=c(2, 3, 4, 1)),
    by=c("time", "age", "sex", "region"),
    dimension_names = dimension_names
)
sexrat_matrix <- afscOM::generate_param_matrix(0.5, dimension_names = dimension_names)

fish_selex_matrix = afscOM::generate_param_matrix(
    aperm(spatial_sablefish_estimates$fish_sel, perm=c(2, 3, 4, 1, 5)),
    by=c("time", "age", "sex", "region", "fleet"),
    dimension_names = dimension_names
)
fish_ret_matrix = afscOM::generate_param_matrix(1, dimension_names = dimension_names, include_fleet_dim = TRUE)
fish_dmr_matrix = afscOM::generate_param_matrix(0, dimension_names = dimension_names, include_fleet_dim = TRUE)
surv_selex_matrix = afscOM::generate_param_matrix(
    aperm(spatial_sablefish_estimates$srv_sel, perm=c(2, 3, 4, 1, 5)),
    by=c("time", "age", "sex", "region", "fleet"),
    dimension_names = c(dimension_names[1:4], list("fleet"=c("LL", "TRWL")))
)

movement <- spatial_sablefish_estimates$Movement
dimnames(movement) <- c(rep(dimension_names[4], 2), dimension_names[1:3])

dem_params <- list(
    waa=waa_matrix,
    mat=mat_matrix,
    mort=mort_matrix,
    sexrat=sexrat_matrix,
    sel=fish_selex_matrix,
    ret=fish_ret_matrix,
    dmr=fish_dmr_matrix,
    surv_sel=surv_selex_matrix,
    movement=movement
)

# dp_y <- subset_dem_params(dem_params, 1, d=1, drop=FALSE)

# 2. Define initial population state
init_naa <- aperm(spatial_sablefish_estimates$NAA[,1,,,drop=FALSE], perm=c(2, 3, 4, 1))
dimnames(init_naa) = c(list("time"=1), dimension_names[2:4])

# 3. Define recruitment timeseries
recruitment <- aperm(spatial_sablefish_estimates$Rec, perm=c(2, 1))
dimnames(recruitment) <- dimension_names[c(1, 4)]

# 4. Define removals timeseries
fmort <- aperm(spatial_sablefish_estimates$Fmort, perm=c(2, 3, 1))
catch <- aperm(spatial_sablefish_estimates$PredCatch, perm=c(2, 3, 1))

# 5. Set Model Options
model_options$removals_input = "F"
model_options$simulate_observations = TRUE
model_options$do_recruits_move = FALSE

model_options$obs_pars <- list(
    # longline fishery, trawl fishery, longline survey, trawl survey
    is_survey   = c(0, 0, 1, 1),  # is this a survey (1) or fishery (0)
    qs          = c(1, 1, 6.113, 9.190), # catchability coefficient (q) for surveys
    catch_cv    = c(0.1, 0.1, 0, 0), # CV on catch for catch observations
    rpn         = c(0, 0, 1, 1), # should RPNs be computed (yes=1, no=0)
    rpn_cv      = c(0, 0, 0.1, 0.1), # RPN CV
    rpw         = c(0, 0, 1, 1), # should RPWs be computed (yes=1, no=0)
    rpw_cv      = c(0, 0, 0.1, 0.1), # RPW CV
    acs         = c(1, 1, 1, 1), # should age compositions be computed (yes=1, no=0)
    ac_samps    = c(200, 100, 200, 100), # total sample size for age composition observations
    ac_as_integers  = c(TRUE, TRUE, TRUE, TRUE), # return age comps as integers (TRUE) or proportions (FALSE)
    acs_agg_sex     = c(FALSE, FALSE, FALSE, FALSE) # should age comps be aggregated by sex
)

om_sim <- afscOM::project(
    init_naa = init_naa, 
    removals_timeseries = fmort, 
    recruitment = recruitment, 
    dem_params = dem_params, 
    nyears = nyears, 
    model_options = model_options
)

ssb <- afscOM::compute_ssb(om_sim$naa, dem_params)
bio <- afscOM::compute_bio(om_sim$naa, dem_params)
catch <- afscOM::compute_total_catch(om_sim$caa)
fleet_catch <- afscOM::compute_fleet_catch(om_sim$caa)
f <- afscOM::compute_total_f(om_sim$faa)

library(patchwork)

# Plot SSB, catch, biomass, and fishing mortality
p1 <- afscOM::plot_ssb(ssb)
p2 <- afscOM::plot_catch(catch)
p3 <- afscOM::plot_bio(bio)
p4 <- afscOM::plot_f(f)
(p1+p2)/(p3+p4)+plot_layout(guides="collect")

matt_ssb <- aperm(spatial_sablefish_estimates$SSB, c(2, 1))
matt_catch <- apply(aperm(spatial_sablefish_estimates$PredCatch, c(2, 3, 1)), c(1, 3), sum)
matt_bio <- apply(
    aperm(spatial_sablefish_estimates$NAA[,1:nyears,,,drop=FALSE], c(2, 3, 4, 1))*dem_params$waa,
    c(1, 4),
    sum
)
matt_f <- apply(aperm(spatial_sablefish_estimates$FAA, c(2, 3, 4, 1, 5)), c(1, 4), max)
dimnames(f) = list("time"=1:nyears, "region"=dimension_names[["region"]])
dimnames(matt_f) = dimnames(f)

afscOM::plot_ssb(ssb, matt_ssb)
afscOM::plot_catch(catch, matt_catch)
afscOM::plot_bio(bio, matt_bio)
afscOM::plot_f(f, matt_f)

# Aggregate survey indices and catches
apply(om_sim$survey_obs$rpns[,,,,1], 1, sum)
apply(om_sim$survey_obs$rpns[,,,,2], 1, sum)
apply(om_sim$land_caa, c(1, 5), sum)

# Aggregate composition data
rpw_props <- t(apply(om_sim$survey_obs$rpws[,,,,2], 1, function(x) x/sum(x)))

apply(om_sim$survey_obs$acs, c(1, 2, 3, 5), sum)

acs <- apply(round(sweep(om_sim$survey_obs$acs, 1, rpw_props, FUN="*"), 0), c(1, 2, 3, 5), sum)
acs

apply(apply(acs, c(1, 4), sum),2, mean)


waa_matrix

nyears <- 200
dimension_names[["time"]] <- 1:nyears

waa_matrix_big <- array(NA, dim=c(nyears, nages, nsexes, nregions))
waa_matrix_big[1:dim(waa_matrix)[1],,,] <- waa_matrix
waa_matrix_big[(dim(waa_matrix)[1]+1):nyears,,,] <- aperm(array(
    rep(waa_matrix[dim(waa_matrix)[1],,,,drop=FALSE], nyears-dim(waa_matrix)[1]),
    dim=c(nages, nsexes, nyears-dim(waa_matrix)[1], nregions)
), c(3, 1, 2, 4))


waa_matrix_big <- extend_array_years(waa_matrix, 200)
dimnames(waa_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:4])
mat_matrix_big <- extend_array_years(mat_matrix, 200)
dimnames(mat_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:4])
mort_matrix_big <- extend_array_years(mort_matrix, 200)
dimnames(mort_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:4])
sexrat_matrix_big <- extend_array_years(sexrat_matrix, 200)
dimnames(sexrat_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:4])
fish_selex_matrix_big <- extend_array_years2(fish_selex_matrix, 200)
dimnames(fish_selex_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:5])
fish_ret_matrix_big <- extend_array_years2(fish_ret_matrix, 200)
dimnames(fish_ret_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:5])
fish_dmr_matrix_big <- extend_array_years2(fish_dmr_matrix, 200)
dimnames(fish_dmr_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:5])
surv_selex_matrix_big <- extend_array_years2(surv_selex_matrix, 200)
dimnames(surv_selex_matrix_big) <- c(list("time"=1:new_years), dimension_names[2:5])
movement_matrix_big <- extend_movement_years(movement, 200)
dimnames(movement_matrix_big) <- c(dimension_names[4], dimension_names[4], list("time"=1:new_years), dimension_names[2], dimension_names[3])


dem_params_big <- list(
    waa = waa_matrix_big,
    mat = mat_matrix_big,
    mort = mort_matrix_big,
    sexrat=sexrat_matrix_big,
    sel = fish_selex_matrix_big,
    ret = fish_ret_matrix_big,
    dmr = fish_dmr_matrix_big,
    surv_sel = surv_selex_matrix_big,
    movement = movement_matrix_big
)

saveRDS(dem_params_big, file=file.path(here::here(), "data", "spatial_dem_params_big.RDS"))

load(file.path(here::here(), "data", "spatial_sablefish_inputs.rda"))
dput(spatial_sablefish_inputs, file.path(here::here(), "data", "spatial_sablefsh_inputs.rdat"))

movement


extend_array_years <- function(og_array, new_years){
    big_array = array(NA, dim=c(new_years, nages, nsexes, nregions))
    big_array[1:dim(og_array)[1],,,] <- og_array
    big_array[(dim(og_array)[1]+1):nyears,,,] <- aperm(array(
        rep(og_array[dim(og_array)[1],,,,drop=FALSE], new_years-dim(og_array)[1]),
        dim=c(nages, nsexes, new_years-dim(og_array)[1], nregions)
    ), c(3, 1, 2, 4))
    # dimnames(big_array) <- c(list("time"=1:new_years), dimension_names[2:4])
    names(dim(big_array)) <- c("nyears", "nages", "nsexes", "nregions")
    return(big_array)
}

extend_array_years2 <- function(og_array, new_years){
    big_array = array(NA, dim=c(new_years, nages, nsexes, nregions, nfleets))
    big_array[1:dim(og_array)[1],,,,] <- og_array
    big_array[(dim(og_array)[1]+1):nyears,,,,] <- aperm(array(
        rep(og_array[dim(og_array)[1],,,,,drop=FALSE], new_years-dim(og_array)[1]),
        dim=c(nages, nsexes, nregions, nfleets, new_years-dim(og_array)[1])
    ), c(5, 1, 2, 3, 4))
    dimnames(big_array) <- c(list("time"=1:new_years), dimension_names[2:5])
    names(dim(big_array)) <- c("nyears", "nages", "nsexes", "nregions", "nfleets")
    return(big_array)
}

extend_movement_years <- function(og_array, new_years){
    big_array = array(NA, dim=c(nregions, nregions, new_years, nages, nsexes))
    big_array[,,1:dim(og_array)[3],,] <- og_array
    big_array[,,(dim(og_array)[3]+1):new_years,,] <- aperm(array(
        rep(og_array[,,dim(og_array)[3],,,drop=FALSE], new_years-dim(og_array)[3]),
        dim=c(nregions, nregions, nages, nsexes, new_years-dim(og_array)[3])
    ), c(1, 2, 5, 3, 4))
    dimnames(big_array) <- c(dimension_names[4], dimension_names[4], list("time"=1:new_years), dimension_names[2], dimension_names[3])
    names(dim(big_array)) <- c("nregions", "nregions", "nyears", "nages", "nsexes")
    return(big_array)
}

