rm(list=ls())

# remotes::install_github('BenWilliams-NOAA/afscOM')
#library(afscOM)
library(devtools)
library(patchwork)
library(tidyverse)
library(ggdist)
library(ggh4x)
library(SpatialSablefishAssessment)

afscOM_dir <- "~/Desktop/Projects/afscOM"
devtools::load_all(afscOM_dir)
source("R/reference_points.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/age_structure_stats.R")
source("R/data_utils.R")
source("R/data_processing.R")
source("R/format_em_data.R")
source("R/fit_TMB_model.R")

assessment <- dget("data/sablefish_assessment_2023.rdat")

# Load OM parameters into global environment
list2env(readRDS("data/sablefish_om.RDS"), globalenv())
model_options$obs_pars$surv_ll$ac_samps <- 50
model_options$obs_pars$surv_ll$rpn_cv <- 0.01
model_options$obs_pars$surv_ll$rpw_cv <- 0.01
model_options$obs_pars$surv_ll$as_integers <- TRUE

model_options$obs_pars$surv_tw$ac_samps <- 50
model_options$obs_pars$surv_tw$rpw_cv <- 0.01
model_options$obs_pars$surv_tw$as_integers <- TRUE

model_options$obs_pars$fish_fx$ac_samps <- 50
model_options$obs_pars$fish_fx$as_integers <- TRUE

model_options$simulate_observations <- TRUE

# Load OM dimensions into global environment
list2env(afscOM::get_model_dimensions(dem_params$sel), globalenv())
nsims <- 1

# Dimension names
dimension_names <- list(
    "time" = 1:nyears,
    "age"  = 2:31,
    "sex"  = c("F", "M"),
    "region" = "alaska",
    "fleet" = c("Fixed", "Trawl"),
    "sims" = 1:nsims
)



# Setup TAC storage
TACs <- rep(0, nyears)
hcr_F <- rep(0, nyears)
out_f <- rep(0, nyears) # vector to store outputted F
TACs[1:64] <- (assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"])


#' 6. Setup empty array to collect derived quantities from the OM
#' It is left to the user to decide what information to store, and
#' in what format they would like to store it.
#'
#' Here, we are storing all OM outputs as they are returned from the
#' OM.
land_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
disc_caa    = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
caa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
faa         = array(NA, dim=c(nyears, nages, nsexes, nregions, nfleets, nsims), dimnames=dimension_names)
tac         = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sims"=1:nsims))
hcr_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sims"=1:nsims))
out_f       = array(NA, dim=c(nyears, 1, 1, 1, nsims), dimnames=list("time"=1:nyears, 1, 1, "region"="Alaska", "sims"=1:nsims))
naa         = array(NA, dim=c(nyears+1, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears+1), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=1:nsims))
naa[1,,,,] = init_naa

naa_est     = array(NA, dim=c(nyears, nages, nsexes, nregions, nsims), dimnames=list("time"=1:(nyears), "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "sim"=1:nsims))

survey_obs <- list(
    ll_rpn = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    tw_rpw = array(NA, dim=c(nyears, 1, 1, nregions)),
    ll_acs = array(NA, dim=c(nyears, nages, nsexes, nregions)),
    fxfish_acs = array(NA, dim=c(nyears, nages, nsexes, nregions))
)

set.seed(1120)
seeds <- sample(1:(nsims*1000), size=nsims, replace=FALSE)
ref_as <- readRDS("data/agestruct_f40.RDS")
ref_naa <- ref_as$naa
for(s in 1:nsims){
    print(s)
    set.seed(seeds[s])
    recruitment <- assessment$natage.female[,1]*2
    projected_recruitment <- sample(recruitment, size=nyears-length(recruitment)+1, replace=TRUE)
    recruitment <- c(recruitment, projected_recruitment)

    for(y in 1:67){
        print(paste("Year", y))
        # Subset the demographic parameters list to only the current year
        # and DO NOT drop lost dimensions.
        dp.y <- subset_dem_params(dem_params = dem_params, y, d=1, drop=FALSE)
        removals_input <- TACs[y]
        fleet.props <- unlist(lapply(model_options$fleet_apportionment, \(x) x[y]))

        prev_naa <- afscOM::subset_matrix(naa[y,,,,s, drop = FALSE], 1, 5, drop=TRUE)
        out_vars <- project(
            removals = removals_input,
            dem_params=dp.y,
            prev_naa=prev_naa,
            recruitment=recruitment[y+1],
            fleet.props = fleet.props,
            options=model_options
        )

        # update state
        land_caa[y,,,,,s] <- out_vars$land_caa_tmp
        disc_caa[y,,,,,s] <- out_vars$disc_caa_tmp
        caa[y,,,,,s] <- out_vars$caa_tmp
        faa[y,,,,,s] <- out_vars$faa_tmp
        naa[y+1,,,,s] <- out_vars$naa_tmp
        out_f[y,,,,s] <- sum(out_vars$F_f_tmp[1,1,1,1,])
        tac[y,,,,s] <- TACs[y]
        hcr_f[y,,,,s] <- hcr_F[y]

        survey_obs$ll_rpn[y,,,] <- out_vars$surv_obs$ll_rpn
        survey_obs$ll_rpw[y,,,] <- out_vars$surv_obs$ll_rpw
        survey_obs$tw_rpw[y,,,] <- out_vars$surv_obs$tw_rpw
        survey_obs$ll_acs[y,,,] <- out_vars$surv_obs$ll_ac_obs
        survey_obs$fxfish_acs[y,,,] <- out_vars$surv_obs$fxfish_caa_obs

        #new_F <- 0.085
        if((y+1) > 64){

            assess_inputs <- format_em_data(
                nyears = y,
                dem_params = dp.y,
                land_caa = out_vars$land_caa_tmp,
                survey_indices = out_vars$surv_obs,
                fxfish_caa_obs = survey_obs$fxfish_acs[y-1,,,,drop=FALSE],
                ll_ac_obs = survey_obs$ll_acs[y-1,,,,drop=FALSE],
                model_options = model_options,
                added_years = 1
            )

            mod_report <- fit_TMB_model(assess_inputs$new_data, assess_inputs$new_parameters)  
            assessment_ssb <- SpatialSablefishAssessment::get_SSB(mod_report) %>% filter(Year == max(Year)) %>% pull(SSB)

            if(y == 64){
                naaf <- t(mod_report$natage_f[,1:63])
                naam <- t(mod_report$natage_m[,1:63])
                naa_est[1:63,,1,,s] <- naaf
                naa_est[1:63,,2,,s] <- naam
            }

            naa_est[y,,1,,s] <- mod_report$natage_f[,y]
            naa_est[y,,2,,s] <- mod_report$natage_m[,y]
            
            ###########

            joint_self <- apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum))
            joint_selm <- apply(dp.y$sel[,,2,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$sel[,,1,,,drop=FALSE], c(1, 2), sum))
            joint_ret <- apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(dp.y$ret[,,1,,,drop=FALSE], c(1, 2), sum))
            
            # reference points are all female based
            ref_pts <- calculate_ref_points(
                nages=nages,
                mort = dp.y$mort[,,1,],
                mat = dp.y$mat[,,1,],
                waa = dp.y$waa[,,1,],
                sel =  joint_self,
                ret = joint_ret,
                avg_rec = mean(recruitment)/2
            )

            #ssb <- apply(out_vars$naa_tmp[,,1,]*dp.y$waa[,,1,,drop=FALSE]*dp.y$mat[,,1,,drop=FALSE], 1, sum)
            hcr_F[y] <- npfmc_tier3_F(assessment_ssb, ref_pts$B40, ref_pts$F40)
            # hcr_F[y] <- as_scalar_threshold_f(
            #                 ssb/ref_pts$B40, 
            #                 naa=out_vars$naa_tmp[,,1,], 
            #                 ref_naa=ref_naa,
            #                 as_func = shannon_diversity,
            #                 #ages = 2:31,
            #                 f_min=0,
            #                 f_max=ref_pts$F40,
            #                 lrp=0.05,
            #                 urp=1.0
            #             )

            joint_sel <- array(NA, dim=dim(out_vars$naa_tmp))
            joint_sel[,,1,] <- joint_self
            joint_sel[,,2,] <- joint_selm

            TACs[y+1] <- simulate_TAC(hcr_F[y], out_vars$naa_tmp, mean(recruitment)/2, joint_sel, dp.y)
        }   
    }

    file.remove("data/sablefish_em_data_curr.RDS")
    file.remove("data/sablefish_em_par_curr.RDS")
}


# l <- afscOM::listN(dem_params, model_options, land_caa, caa, naa, faa, tac, hcr_f, survey_obs)
# saveRDS(l, "data/om2.RDS")

library(reshape2)
library(ggdist)
library(ggh4x)

# sum(naa_est[70, , 1, 1,1]*dp.y$mat[,,1,1]*dp.y$waa[,,1,1])
# sum(naa[70, , 1, 1,1])

om_ssb <- apply(naa[1:100,,1,1,1]*dem_params$mat[1:100,,1,1]*dem_params$waa[1:100,,1,1], 1, sum)
em_ssb <- apply(naa_est[1:100,,1,1,1]*dem_params$mat[1:100,,1,1]*dem_params$waa[1:100,,1,1], 1, sum)

plot(1:100, om_ssb, type="l", ylim=c(0, 300))
lines(1:100, em_ssb, col="red")

pos_scales <- list(
    scale_y_continuous(limits=c(0, 300)),
    scale_y_continuous(limits=c(0, 100)),
    scale_y_continuous(limits=c(0, 750)),
    scale_y_continuous(limits=c(0, 0.2))
)

biomass_df <- create_summary_biomass_df(naa, dem_params$waa, dem_params$mat)

catch_df <- create_summary_catch_df(caa, out_f)

d <- create_biomass_catch_summary_df(biomass_df,catch_df)

ggplot(d)+
    geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper), size=0.5)+
    geom_vline(xintercept=2023-1960+1, linetype="dashed")+
    geom_hline(data=ref_points_df, aes(yintercept=rp), linetype="dashed", color="red")+
    scale_fill_brewer(palette = "Blues")+
    scale_y_continuous(limits=c(0, 100))+
    facet_wrap(~name, scales="free_y")+
    facetted_pos_scales(y=pos_scales)+
    scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
    coord_cartesian(expand=0)+
    theme_bw()+
    theme(
        strip.background = element_blank(),
        strip.text = element_text(size=14, hjust=0),
        panel.spacing.x = unit(0.75, "cm"),
        axis.title = element_blank(),
        legend.position = "bottom"
    )

# ggsave("~/Desktop/tier3_hcr.jpeg")

final_ref_points_tot <- calculate_ref_points(nages, dp.y$mort, dp.y$mat, dp.y$waa, joint_sel, array(1, dim=dim(joint_sel)), mean(recruitment))
final_ref_points_fem <- calculate_ref_points(nages, dp.y$mort[,,1,], dp.y$mat[,,1,], dp.y$waa[,,1,], joint_sel[,,1,], joint_ret, mean(recruitment)/2)

ref_points_df <- data.frame(
    name = factor(c("tot_spbio", "catch", "tot_bio", "F"), levels=c("tot_spbio", "catch", "tot_bio", "F"), labels=c("Spawning Biomass (mt)", "Catch (mt)", "Total Biomass (mt)", "Summary F")),
    rp = c(final_ref_points_fem$B40, final_ref_points_tot$B40*final_ref_points_tot$F40, final_ref_points_tot$B40, final_ref_points_fem$F40)
)


melt(naa) %>% as_tibble() %>%
    filter(sex == "F") %>%
    drop_na() %>%
    group_by(time, age, sim) %>%
    summarise(
        total_num = sum(value)
    ) %>%
    group_by(time, age) %>%
    summarise(
        total_num = mean(total_num)
    ) %>%
    mutate(
        class = case_when(age < 7 ~ "Young", age < 16 ~ "Immature", age > 15 ~ "Mature")
    ) %>%
    mutate(
        class = factor(class, levels=c("Young", "Immature", "Mature"), labels=c("Young (<7yo)", "Immature (7-15yo)", "Mature (16+yo)"))
    ) %>%
    
    ggplot()+
        geom_bar(aes(x=time, y=total_num, fill=class), position="fill", stat="identity")+
        scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
        coord_cartesian(expand=0)
