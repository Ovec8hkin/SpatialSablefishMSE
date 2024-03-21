rm(list=ls())

# remotes::install_github('BenWilliams-NOAA/afscOM')
#library(afscOM)
library(devtools)
library(patchwork)
library(tidyverse)

afscOM_dir <- "~/Desktop/Projects/afscOM"
devtools::load_all(afscOM_dir)
source("R/reference_points.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")

assessment <- dget("data/sablefish_assessment_2023.rdat")

# Load OM parameters into global environment
list2env(readRDS("data/sablefish_om.RDS"), globalenv())

# Load OM dimensions into global environment
list2env(afscOM::get_model_dimensions(dem_params$sel), globalenv())
nsims <- 50

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

seeds <- sample(1:(nsims*1000), size=nsims, replace=FALSE)
for(s in 1:50){
    print(s)
    set.seed(seeds[s])
    recruitment <- assessment$natage.female[,1]*2
    projected_recruitment <- sample(recruitment, size=nyears-length(recruitment)+1, replace=TRUE)
    recruitment <- c(recruitment, projected_recruitment)

    for(y in 1:(nyears-1)){

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

        #new_F <- 0.085
        if((y+1) > 64){
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

            ssb <- apply(out_vars$naa_tmp[,,1,]*dp.y$waa[,,1,,drop=FALSE]*dp.y$mat[,,1,,drop=FALSE], 1, sum)
            hcr_F[y] <- npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)

            joint_sel <- array(NA, dim=dim(out_vars$naa_tmp))
            joint_sel[,,1,] <- joint_self
            joint_sel[,,2,] <- joint_selm

            TACs[y+1] <- simulate_TAC(hcr_F[y], out_vars$naa_tmp, mean(recruitment)/2, joint_sel, dp.y)
        }   
    }
}

library(reshape2)
library(ggdist)
library(ggh4x)


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

ggsave("~/Desktop/tier3_hcr.jpeg")

final_ref_points_tot <- calculate_ref_points(nages, dp.y$mort, dp.y$mat, dp.y$waa, joint_sel, array(1, dim=dim(joint_sel)), mean(recruitment))
final_ref_points_fem <- calculate_ref_points(nages, dp.y$mort[,,1,], dp.y$mat[,,1,], dp.y$waa[,,1,], joint_sel[,,1,], joint_ret, mean(recruitment)/2)

ref_points_df <- data.frame(
    name = factor(c("tot_spbio", "catch", "tot_bio", "F"), levels=c("tot_spbio", "catch", "tot_bio", "F"), labels=c("Spawning Biomass (mt)", "Catch (mt)", "Total Biomass (mt)", "Summary F")),
    rp = c(final_ref_points_fem$B40, final_ref_points_tot$B40*final_ref_points_tot$F40, final_ref_points_tot$B40, final_ref_points_fem$F40)
)
