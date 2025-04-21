rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)
library(SPoCK)

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

devtools::load_all(afscOM_dir)

lapply(list.files("R", full.names = TRUE, include.dirs = FALSE, recursive = TRUE), source)

## Helpful functions -----------------
rss <- function(obs, exp){
    num <- sum(obs*(1-obs))
    denom <- sum((exp-obs)^2)
    return(num/denom)
}

aggregate_comps <- function(y, nfleets, weight_type=1){
    tmp <- array(NA, dim=c(1, 30, 2, 1, nfleets))
    for(f in 1:nfleets){
        ISS <- om$model_options$obs_pars$ac_samps[f]
        agg_sex <- om$model_options$obs_pars$acs_agg_sex[f]
        as_int <- om$model_options$obs_pars$ac_as_integers[f]

        is_survey <- om$model_options$obs_pars$is_survey[f]
        if(is_survey){
            selex <- subset_matrix(dp_y$surv_sel, r=f-2, d=5, drop=TRUE)
            weights <- apply(model_runs$naa[y,,,,drop=FALSE]*om$dem_params$waa[y,,,,drop=FALSE], 4, sum)
        }else{
            selex <- subset_matrix(dp_y$sel, r=f, d=5, drop=TRUE)
            weights <- apply(model_runs$caa[y,,,,,drop=FALSE], 4, sum)
        }

        test <- generate_aggregated_comp(
            ac = model_runs$naa[y,,,,drop=FALSE],
            weight_type = weight_type,
            weights = rep(0.2, 5),
            selex = selex,
            total_samples = ISS,
            aggregate_sex = agg_sex,
            as_integers = as_int
        )

        tmp[,,,,f] <- test
    } 
    return(tmp)
}

generate_aggregated_comp <- function(ac, weight_type, selex, weights, total_samples, aggregate_sex, as_integers){
    if(weight_type == 1){
        agg_comp <- simulate_weighted_comps_ISS(ac, selex, weights, total_samples, aggregate_sex, as_integers)
    }else if(weight_type == 2){
        agg_comp <- simulate_weighted_comps_SAMPLE(ac, selex, weights, total_samples, aggregate_sex, as_integers)
    }else{
        agg_comp <- simulate_weighted_comps_HYRBID(ac, selex, weights, total_samples, aggregate_sex, as_integers)
    }
    return(agg_comp)
}

om_to_spock <- function(x){
    return(aperm(x, perm=c(4, 1, 2, 3)))
}

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons
nyears <- 100
sable_om <- readRDS("data/spatial_sablefish_om.RDS") # Read this saved OM from a file
source("dev/oms.R")
source("dev/hcrs.R")


#' 3. Run the closed-loop MSE simulation
#' A single MSE simulation can be run using the `run_mse(...)`
#' function, while multiple MSE simulations can be run (serially)
#' using the `run_mse_multiple(...)` function.
#' 

mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 62
mse_options$recruitment_start_year <- 62
mse_options$n_proj_years <- 20
mse_options$run_estimation <- FALSE

om <- om_rand_recruit
om$model_options$obs_pars$ac_samps <- c(10, 5, 10, 5)*20
mp <- mp_f40

nsims <- 50
seeds <- sample(1:1e6, nsims)

nyears <- mse_options$n_proj_years+mse_options$n_spinup_years
om_ssbs <- array(NA, dim=c(nyears, nsims))
em_ssbs <- array(NA, dim=c(nyears, nsims))

for(s in 1:nsims){

    model_runs <- run_mse(
        om=om,
        mp=mp,
        mse_options=mse_options,
        seed=seeds[s]
    )

    rISS <- get_rISS(nyears, nsamples=10, weight_type = 3, om, model_runs)

    model_dimensions <- afscOM::get_model_dimensions(om$dem_params$sel)

    aggregated_survey_obs <- aggregate_observations(model_runs$survey_obs, nyears)
    aggregated_acs <- array(NA, dim=c(nyears, 30, 2, 1, 4))
    for(i in 1:nyears){
        aggregated_acs[i,,,,] <- aggregate_comps(i, nfleets=4, weight_type = 1)
    }
    aggregated_survey_obs$acs <- aggregated_acs

    input_list <- generate_RTMB_inputs(nyears, om, model_runs, aggregated_survey_obs, rISS)

    data <- input_list$data
    parameters <- input_list$par
    mapping <- input_list$map

    sabie_rtmb_model <- fit_model(
        data,
        parameters,
        mapping,
        random = NULL,
        newton_loops = 3,
        silent = FALSE
    )

    om_ssb <- apply(model_runs$naa[1:nyears,,1,,drop=FALSE]*om$dem_params$mat[1:nyears,,1,,drop=FALSE]*om$dem_params$waa[1:nyears,,1,,drop=FALSE], 1, sum)
    om_ssbs[,s] <- om_ssb
    em_ssbs[,s] <- sabie_rtmb_model$rep$SSB[1,]

}

sd <- RTMB::sdreport(sabie_rtmb_model)

jit <- do_jitter(
    data,
    parameters,
    mapping,
    sd=0.1,
    n_jitter=50,
    n_newton_loops = 3,
    do_par=TRUE,
    n_cores=15
)

prop_converged <- jit %>% filter(Year == 1, Type == 'SSB') %>% summarize(prop_conv = sum(Hessian) / length(Hessian))
final_mod <- reshape2::melt(sabie_rtmb_model$rep$SSB) %>% 
    rename(Region = Var1, Year = Var2) %>% 
    mutate(Type = 'SSB') %>% 
    bind_rows(
        reshape2::melt(sabie_rtmb_model$rep$Rec) %>% 
        rename(Region = Var1, Year = Var2) %>% 
        mutate(Type = 'Recruitment')
    ) 
 
ggplot() +
     geom_line(jit, mapping = aes(x = Year + 1959, y = value, group = jitter, color = Hessian), lwd = 1) +
     geom_line(final_mod, mapping = aes(x = Year + 1959, y = value), color = "black", lwd = 1.3 , lty = 2) +
     facet_grid(Type~Region, scales = 'free') +
     labs(x = "Year", y = "Value") +
     theme_bw(base_size = 20) +
     scale_color_manual(values = c("red", 'grey')) +
     geom_text(data = jit %>% filter(Type == 'SSB', Year == 1, jitter == 1), 
               aes(x = Inf, y = Inf, label = paste("Proportion Converged: ", round(prop_converged$prop_conv, 3))),
               hjust = 1.1, vjust = 1.9, size = 6, color = "black")

get_rISS <- function(nyears, nsamples, weight_type, om, model_runs){
    age_comp_obs <- array(NA, dim=c(nyears, 30, 2, 1, 4, nsamples))
    r_sampsize <- array(NA, dim=c(nyears, 4, nsamples))
    for(s in 1:nsamples){
        print(s)
        for(y in 1:nyears){
            age_comp_obs[y,,,,,s] <- aggregate_comps(y, nfleets=4, weight_type = weight_type)
            for(f in 1:4){
                if(f < 3){
                    selected <- model_runs$naa[y,,,,drop=FALSE]*subset_matrix(om$dem_params$sel[y,,,,,drop=FALSE], f, 5, drop=TRUE)
                }else{
                    selected <- model_runs$naa[y,,,,drop=FALSE]*subset_matrix(om$dem_params$surv_sel[y,,,,,drop=FALSE], f-2, 5, drop=TRUE)
                }
                
                total_selected <- array(apply(selected, c(1, 2, 3), sum), dim=c(1,30,2,1))
                selected_prop <- aperm(apply(total_selected, c(1), \(x) x/sum(x)), c(2, 1))

                ac <- array(age_comp_obs[y,,1:2,,1,s], c(30*2, 1))
                ac <- ac/sum(ac)

                r_sampsize[y,f,s] <- rss(obs=selected_prop[1,], exp=ac)
            }
        }
    }

    r_ISS <- round(apply(r_sampsize, c(1, 2), mean), 0)
    return(r_ISS)
}



bind_rows(process_ssb_trajectories(listN(om_ssbs, em_ssbs)) %>% mutate(om="Age Movement", weight="ISS")) %>%
    ggplot()+
        geom_lineribbon(aes(x=Year, y=EM_SSB, ymin=EM_SSB.lower, ymax=EM_SSB.upper), linewidth=0.1)+
        geom_pointinterval(aes(x=Year, y=OM_SSB, ymin=OM_SSB.lower, ymax=OM_SSB.upper), color="red", size=1)+
        scale_fill_brewer(palette="Blues")+
        coord_cartesian(ylim=c(0, 300))+
        facet_grid(rows=vars(weight), cols=vars(om))+
        custom_theme


om_ssb <- apply(model_runs$naa[1:nyears,,1,,drop=FALSE]*om$dem_params$mat[1:nyears,,1,,drop=FALSE]*om$dem_params$waa[1:nyears,,1,,drop=FALSE], 1, sum)
plot(om_ssb,type="l", ylim=c(0, 300))
lines(sabie_rtmb_model$rep$SSB[1,], col="red")
