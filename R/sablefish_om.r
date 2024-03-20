rm(list=ls())

# remotes::install_github('BenWilliams-NOAA/afscOM')
#library(afscOM)
library(devtools)
library(patchwork)
library(tidyverse)

devtools::load_all("~/Desktop/Side projects/afscOM")
source("~/Desktop/Side projects/SablefishMSE/R/reference_points.R")
source("~/Desktop/Side projects/SablefishMSE/R/harvest_control_rules.R")

assessment <- dget("~/Desktop/Side projects/afscOM/data/test.rdat")

#' 1. Define model dimensions and dimension names
#'  - nyears: number of years over which to simulate
#'  - nages: number of age classes to use
#'  - nsexes: number of sexes (1 = females only)
#'  - nregions: number of spatial regions
#'  - nfleets: number of fishing fleets (survey fleets are handled seperately)
#'
#' Dimension names must be defined in order to use the helper function
#' `generate_param_matrix` which handles filling complex multi-dimensional
#' arrays with smaller dimensions values, vectors, or matrices.
nyears <- 101
nages  <- 30
nsexes <- 2
nregions <- 1
nfleets <- 2
nsurveys <- 2

dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 2:31,
    "sex"  = c("F", "M"),
    "region" = "alaska",
    "fleet" = c("Fixed", "Trawl")
)

model_params <- set_model_params(nyears, nages, nsexes, nregions, nfleets)

#' 2. Generate list with appropriate demographic parameters
#' Required parameters are:
#'  - waa: weight-at-age
#'  - mat: maturity-at-age
#'  - mort: natural-mortality-at-age (instantaneous)
#'  - sexrat: sex ratio
#'  - sel: selectivity-at-age
#'  - ret: retention-at-age
#'  - dmr: discard-mortality-at-age (instantaneous)
#'  - surv_sel: survey selectivity-at-age

M <- 0.113179
mort <- generate_param_matrix(M, dimension_names = dimension_names)
mort[,,2,] <- mort[,,2,]-0.00819813

prop_males <- 0.5
sexrat <- generate_param_matrix(prop_males, dimension_names = dimension_names)

retention <- 1.0
ret <- generate_param_matrix(retention, dimension_names = dimension_names, include_fleet_dim = TRUE)

discard <- 0.0
dmr <- generate_param_matrix(discard, dimension_names = dimension_names, include_fleet_dim = TRUE)

maturity <- assessment$growthmat[, "mage.block1"]
mat <- generate_param_matrix(maturity, dimension_names = dimension_names, by="age")

weight_mat <- matrix(NA, nrow=nages, ncol=nsexes, dimnames=dimension_names[c("age", "sex")])
weight_mat[,1] <- assessment$growthmat[, "wt.f.block1"]
weight_mat[,2] <- assessment$growthmat[, "wt.m.block1"]
waa <- generate_param_matrix(weight_mat, dimension_names = dimension_names, by=c("age", "sex"))

selex_mat <- array(NA, dim=c(nages, nsexes, nfleets), dimnames=dimension_names[c("age", "sex", "fleet")])
selex_mat[,1,1] <- assessment$agesel[,"fish1sel.f"]
selex_mat[,2,1] <- assessment$agesel[,"fish1sel.m"]
selex_mat[,1,2] <- assessment$agesel[,"fish3sel.f"]
selex_mat[,2,2] <- assessment$agesel[,"fish3sel.m"]
sel <- generate_param_matrix(selex_mat, dimension_names = dimension_names, by=c("age", "sex", "fleet"), include_fleet_dim = TRUE)
sel[(36:56),,1,,1] <- matrix(rep(assessment$agesel[, "fish4sel.f"], length(36:56)), ncol=nages, byrow=TRUE)
sel[(36:56),,2,,1] <- matrix(rep(assessment$agesel[, "fish4sel.m"], length(36:56)), ncol=nages, byrow=TRUE)
sel[(57:nyears),,1,,1] <- matrix(rep(assessment$agesel[, "fish5sel.f"], length(57:nyears)), ncol=nages, byrow=TRUE)
sel[(57:nyears),,2,,1] <- matrix(rep(assessment$agesel[, "fish5sel.m"], length(57:nyears)), ncol=nages, byrow=TRUE)

survey_selex_mat <- array(NA, dim=c(nages, nsexes, 2), dimnames=dimension_names[c("age", "sex", "fleet")])
survey_selex_mat[,1,1] <- assessment$agesel[,"srv1sel.f"]
survey_selex_mat[,2,1] <- assessment$agesel[,"srv1sel.m"]
survey_selex_mat[,1,2] <- assessment$agesel[,"srv7sel.f"]
survey_selex_mat[,2,2] <- assessment$agesel[,"srv7sel.m"]
survey_sel <- generate_param_matrix(survey_selex_mat, dimension_names = dimension_names, by=c("age", "sex", "fleet"), include_fleet_dim = TRUE)
survey_sel[(57:nyears),,1,,1] <- matrix(rep(assessment$agesel[, "srv10sel.f"], length(57:nyears)), ncol=nages, byrow=TRUE)
survey_sel[(57:nyears),,2,,1] <- matrix(rep(assessment$agesel[, "srv10sel.m"], length(57:nyears)), ncol=nages, byrow=TRUE)

dem_params <- list(
    waa=waa,
    mat=mat,
    mort=mort,
    sexrat=sexrat,
    sel=sel,
    ret=ret,
    dmr=dmr,
    surv_sel=survey_sel
)

#dem_params <- validate_dem_params(dem_params, model_params)

#' 3. Define starting population condition (NAA)
#' Here, we will use the NAA in 1960 from the 2023
#' Alaska sablefish assessment.
init_naa <- array(NA, dim=c(1, nages, nsexes, nregions), dimnames = list(1, 2:31, c("F", "M"), "alaska"))
init_naa[,,1,] <- assessment$natage.female["1960",]
init_naa[,,2,] <- assessment$natage.male["1960",]

#' 4. Define recruitment timeseries
#' Where possible (e.g. where recruitment is independent of
#' the population state), recruitment should be defined
#' external to the OM simulation loop to ensure reproducability
#' across model runs.
#'
#' Here we using the recruitment timeseries from the 2023 Alaska
#' sablefish assessment.

# Multiply by 2 because two sexes
recruitment <- assessment$natage.female[,1]*2
projected_recruitment <- sample(recruitment, size=nyears-length(recruitment)+1, replace=TRUE)
recruitment <- c(recruitment, projected_recruitment)


#' 5. Define catch history
#' Catch histories are defined by three variables:
#'  - TAC: the total allowable catch each year across the entire model
#'  - region apportionment: the proportion of the TAC allocated to each
#'                          region in the model
#'  - fleet apportionment: the proportion of the TAC allocation to each
#'                         fishing fleet in each region in the model
#' The apportionment timeseries are placed in a `model_options` list.
#'
#' Alternatively, a timeseries of fishing mortality rates, F, can be
#' provided in lieu of TACs. In this case, the timeseries of F for each
#' fleet, in each area, must be provided in an array of dimensions
#' [nyears, 1, 1, nregions, nfleets].
#'
#' Always specify the `removals_input` in the `model_options` list.
#' `removals_input` can be either "catch" or "F", depending on what
#' data you are providing as removals input.
TACs <- rep(0, nyears)
hcr_F <- rep(0, nyears)
out_f <- rep(0, nyears) # vector to store outputted F
TACs[1:64] <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])

fixed_fleet_prop <- rep(NA, nyears)
fixed_fleet_prop[1:64] <- assessment$t.series[,"Catch_HAL"]/(assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])
fixed_fleet_prop[65:nyears] <- mean(fixed_fleet_prop[1:64])
trawl_fleet_prop <- 1-fixed_fleet_prop

model_options <- list(
    region_apportionment = list(1),
    fleet_apportionment = list(fixed_fleet_prop, trawl_fleet_prop),
    removals_input = "catch"
)


# f_timeseries <- assessment$t.series[,c("F_HAL", "F_TWL")] %>% as.matrix
# f_timeseries <- array(f_timeseries, dim=c(nyears, 1, 1, 1, 2),
#                 dimnames = list("time"=1:nyears,
#                                 age="all",
#                                 sex="all",
#                                 "region"="alaska",
#                                 "fleet"=c("Fixed", "Trawl")))

# model_options <- list(
#     removals_input = "F"
# )

#' 6. Define parameters for observation processes
#' Observation process parameters include catchability coefficients
#' (q), observation errors, and sample sizes for age/length comps.

obs_pars <- list(
    surv_ll = list(
        q = 6.41338,
        rpn_cv = 0.20,
        rpw_cv = 0.10,
        ac_samps = 1000
    ),
    surv_tw = list(
        q = 0.8580,
        rpw_cv = 0.10,
        ac_samps = 1000
    ),
    fish_fx = list(
        ac_samps = 1000
    )
)

model_options$obs_pars <- obs_pars

#' 6. Setup empty array to collect derived quantities from the OM
#' It is left to the user to decide what information to store, and
#' in what format they would like to store it.
#'
#' Here, we are storing all OM outputs as they are returned from the
#' OM.
land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets))
naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions))
naa[1,,,] = init_naa

survey_preds <- list(
    ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_ac = array(NA, dim=c(nyears, nages, nsexes, nregions)),
    fxfish_caa = array(NA, dim=c(nyears, nages, nsexes, nregions))
)

survey_obs <- list(
    ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions)),
    fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions))
)

#' 7. Run the OM forward in time
#' The `project` function handles projecting the population
#' forward in time 1 timestep based on the provded TAC, current
#' population size, provided recruitment, and other model
#' options. In order to simulate multiple years, the `project`
#' function MUST be wrapped in a loop over the intended number
#' of simulation years.
#'
#' The `project` function does not "know" about the presence of
#' an external "year loop", and thus expects arguments to
#' represent data for a single year at a time. Helper functions
#' have been provided to help with correctly subsetting the
#' demographic matrices to ensure input data is of the correct
#' dimensionality.

set.seed(1007)
for(y in 1:(nyears-1)){

    # Subset the demographic parameters list to only the current year
    # and DO NOT drop lost dimensions.
    dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
    removals_input <- TACs[y]
    fleet.props <- unlist(lapply(model_options$fleet_apportionment, \(x) x[y]))
    out_vars <- project(
        removals = removals_input,
        dem_params=dp.y,
        prev_naa=naa[y,,,, drop = FALSE],
        recruitment=recruitment[y+1],
        fleet.props = fleet.props,
        options=model_options
    )

    # update state
    land_caa[y,,,,] <- out_vars$land_caa_tmp
    disc_caa[y,,,,] <- out_vars$disc_caa_tmp
    caa[y,,,,] <- out_vars$caa_tmp
    faa[y,,,,] <- out_vars$faa_tmp
    naa[y+1,,,] <- out_vars$naa_tmp
    out_f[y] <- sum(out_vars$F_f_tmp[1,1,1,1,])

    survey_preds$ll_rpn[y,,,] <- out_vars$surv_preds$ll_rpn
    survey_preds$ll_rpw[y,,,] <- out_vars$surv_preds$ll_rpw
    survey_preds$tw_rpw[y,,,] <- out_vars$surv_preds$tw_rpw
    survey_preds$ll_ac[y,,,] <- out_vars$surv_preds$ll_ac
    survey_preds$fxfish_caa[y,,,] <- out_vars$surv_preds$fxfish_caa

    survey_obs$ll_rpn[y,,,] <- out_vars$surv_obs$ll_rpn
    survey_obs$ll_rpw[y,,,] <- out_vars$surv_obs$ll_rpw
    survey_obs$tw_rpw[y,,,] <- out_vars$surv_obs$tw_rpw
    survey_obs$ll_acs[y,,,] <- out_vars$surv_obs$ll_ac_obs
    survey_obs$fxfish_acs[y,,,] <- out_vars$surv_obs$fxfish_caa_obs

    #new_F <- 0.085
    if((y+1) > 64){
        joint_self <- apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum))
        joint_selm <- apply(dp.y$sel[,,2,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum))
        joint_ret <- apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum))
        F35 <- spr_x(
            nages=30,
            mort = dp.y$mort[,,1,],
            mat = dp.y$mat[,,1,],
            waa = dp.y$waa[,,1,],
            sel =  joint_self,
            ret = joint_ret,
            target_x = 0.35
        )
        F40 <- spr_x(
            nages=30,
            mort = dp.y$mort[,,1,],
            mat = dp.y$mat[,,1,],
            waa = dp.y$waa[,,1,],
            sel =  joint_self,
            ret = joint_ret,
            target_x = 0.40
        )
        B40 <- compute_bx(
            nages=30,
            mort = dp.y$mort[,,1,],
            mat = dp.y$mat[,,1,],
            waa = dp.y$waa[,,1,],
            sel =  joint_self,
            ret = joint_ret,
            F = F40,
            avg_rec = mean(recruitment/2)
        )

        ssb <- apply(out_vars$naa_tmp[,,1,]*dp.y$waa[,,1,,drop=FALSE]*dp.y$mat[,,1,,drop=FALSE], 1, sum)
        hcr_F[y] <- npfmc_tier3_F(ssb, B40, F40, 0.05)

        # calculate quantities to get TAC
        # Project population forward by 1 year with terminal F
        proj_N <- array(NA, dim = c(dim(out_vars$naa_tmp)))
        jointsel <- rbind(joint_self, joint_selm)

        for(s in 1:nsexes) {
          for(a in 1:(nages - 1)) {
            if(a == 1) proj_N[,1,s,] <- mean(recruitment) * 0.5 # using this for now, but assessment uses 1979 - terminal
            else proj_N[,a,s,] <- out_vars$naa_tmp[,a,s,] * exp((out_f[y] * jointsel[s,a]) + M)
          } # end a
          proj_N[1,nages,s,1] <- out_vars$naa_tmp[1,nages-1,s,1]*exp(-(out_f[y] * jointsel[s,nages-1]) + M) +
                                 out_vars$naa_tmp[1,nages,s,1]*exp(-(out_f[y] * jointsel[s,nages]) + M)
          # Now calculate TACs
          TAC_tmp <- 0
          FAA_tmp <- (hcr_F[y] * jointsel[s,])
          ZAA_tmp <- FAA_tmp + M
          TAC_tmp <- sum(FAA_tmp / ZAA_tmp * proj_N[,,s,] * (1 - exp(-ZAA_tmp)) * dp.y$waa[,,s,])
          TAC_tmp <- TAC_tmp + TAC_tmp # update TAC for females and male
        } # end s

        TACs[y+1] <- TAC_tmp
    }

}


#' 8. Plot OM Results
#'

ssb <- apply(naa[1:nyears,,1,]*dem_params$waa[,,1,]*dem_params$mat[,,1,], 1, sum)
bio <- apply(naa[1:nyears,,,]*dem_params$waa[,,,], 1, sum)
catch <- apply(caa, 1, sum)
# f <- apply(apply(faa, c(1, 5), \(x) max(x)), 1, sum)

ssb_comp <- data.frame(
    year=1960:(1960+nyears-1),
    #assess_ssb=assessment$t.series[, "spbiom"],
    om_ssb=ssb,
    assess_catch=TACs,
    om_catch=catch,
    #assess_bio=assessment$t.series[, "totbiom"],
    om_bio = bio,
    #assess_f = assessment$t.series[,"fmort"],
    om_f = out_f
)

library(ggplot2)

p1 <- ggplot(ssb_comp, aes(x=year))+
    #geom_line(aes(y=assess_ssb, color="Assessment"))+
    geom_line(aes(y=om_ssb, color="OM"), size=0.7)+
    geom_hline(yintercept = B40) +
    scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="SSB", x="Year", title="Spawning Biomass Comparison")+
    theme_bw()

p2 <- ggplot(ssb_comp, aes(x=year))+
    #geom_line(aes(y=assess_catch, color="Assessment"))+
    geom_line(aes(y=om_catch, color="OM"), size=0.7)+
    scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 10))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="Catch", x="Year", title="Total Catch Comparison")+
    theme_bw()

p3 <- ggplot(ssb_comp, aes(x=year))+
    #geom_line(aes(y=assess_bio, color="Assessment"))+
    geom_line(aes(y=om_bio, color="OM"), size=0.7)+
    scale_y_continuous(limits=c(0, 1000), breaks=seq(0, 1000, 100))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="Biomass", x="Year", title="Total Biomass Comparison")+
    theme_bw()

p4 <- ggplot(ssb_comp %>% filter(om_f != 0), aes(x=year))+
    #geom_line(aes(y=assess_f, color="Assessment"))+
    geom_line(aes(y=om_f, color="OM"), size=0.7)+
    geom_hline(yintercept = F40) +
    scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.2, 0.05))+
    scale_x_continuous(breaks=seq(1960, 2020, 10))+
    coord_cartesian(expand=0)+
    scale_color_manual(name="Model", values=c("black", "red"))+
    labs(y="F", x="Year", title="Fishing Mortality Comparison")+
    theme_bw()

p <- (p1+p2)/(p3+p4)+plot_layout(guides="collect")
p
#ggsave("~/Desktop/sablefish_assess_om.png", plot=p, width=8, height=8, units=c("in"))


ll_surv_data <- data.frame(assessment$obssrv3) %>% rownames_to_column("Year")
ll_surv_data$Year <- as.numeric(ll_surv_data$Year)
ll_surv_data$om_pred <- survey_preds$ll_rpn[31:64,,,]
ll_surv_data$om <- survey_obs$ll_rpn[31:64,,,]
ll_surv_data$om.lci <- ll_surv_data$om -1.96*0.20*ll_surv_data$om
ll_surv_data$om.uci <- ll_surv_data$om +1.96*0.20*ll_surv_data$om

p1 <- ggplot(ll_surv_data, aes(x=Year, y=obssrv3, group=1))+
    geom_pointrange(aes(ymin=obssrv3.lci, ymax=obssrv3.uci, color="Assessment"))+
    geom_line(aes(y=predsrv3, color="Assessment"))+
    geom_line(aes(y=om_pred, color="OM"))+
    geom_pointrange(aes(y=om, ymin=om.lci, ymax=om.uci, color="OM"))+
    scale_y_continuous(limits=c(0, 3000), breaks=seq(0, 3000, 500))+
    scale_x_continuous(limits=c(1960, 2025), breaks=seq(1960, 2023, 10))+
    scale_color_manual(name="Model", values=c("black", "blue"))+
    coord_cartesian(expand=0)+
    labs(y="LL Survey RPN", x="Year", title="LL Survey RPN Comparison")+
    theme_bw()

ll_surv_data <- data.frame(assessment$obssrv1) %>% rownames_to_column("Year")
ll_surv_data$Year <- as.numeric(ll_surv_data$Year)
ll_surv_data$om_pred <- survey_preds$ll_rpw[31:64,,,]
ll_surv_data$om <- survey_obs$ll_rpw[31:64,,,]
ll_surv_data$om.lci <- ll_surv_data$om -1.96*0.10*ll_surv_data$om
ll_surv_data$om.uci <- ll_surv_data$om +1.96*0.10*ll_surv_data$om

p2 <- ggplot(ll_surv_data, aes(x=Year, y=obssrv1, group=1))+
    geom_pointrange(aes(ymin=obssrv1.lci, ymax=obssrv1.uci, color="Assessment"))+
    geom_line(aes(y=predsrv1, color="Assessment"))+
    geom_line(aes(y=om_pred, color="OM"))+
    geom_pointrange(aes(y=om, ymin=om.lci, ymax=om.uci, color="OM"))+
    scale_y_continuous(limits=c(0, 5500), breaks=seq(0, 5500, 1000))+
    scale_x_continuous(limits=c(1960, 2025), breaks=seq(1960, 2023, 10))+
    scale_color_manual(name="Model", values=c("black", "blue"))+
    coord_cartesian(expand=0)+
    labs(y="LL Survey RPW", x="Year", title="LL Survey RPW Comparison")+
    theme_bw()

tw_surv_data <- data.frame(assessment$obssrv7) %>% rownames_to_column("Year")
tw_surv_data$Year <- as.numeric(tw_surv_data$Year)
tw_surv_data$om_pred <- survey_preds$tw_rpw[tw_surv_data$Year-1960+1,,,]
tw_surv_data$om <- survey_obs$tw_rpw[tw_surv_data$Year-1960+1,,,]
tw_surv_data$om.lci <- tw_surv_data$om -1.96*0.10*tw_surv_data$om
tw_surv_data$om.uci <- tw_surv_data$om +1.96*0.10*tw_surv_data$om

p3 <- ggplot(tw_surv_data, aes(x=Year, y=obssrv7, group=1))+
    geom_pointrange(aes(ymin=obssrv7.lci, ymax=obssrv7.uci, color="Assessment"))+
    geom_line(aes(y=predsrv7, color="Assessment"))+
    geom_line(aes(y=om_pred, color="OM"))+
    geom_pointrange(aes(y=om, ymin=om.lci, ymax=om.uci, color="OM"))+
    scale_y_continuous(limits=c(0, 500), breaks=seq(0, 500, 100))+
    scale_x_continuous(limits=c(1960, 2025), breaks=seq(1960, 2023, 10))+
    scale_color_manual(name="Model", values=c("black", "blue"))+
    coord_cartesian(expand=0)+
    labs(y="TW Survey RPW", x="Year", title="TW Survey RPW Comparison")+
    theme_bw()

p1/p2/p3 + plot_layout(guides="collect")

