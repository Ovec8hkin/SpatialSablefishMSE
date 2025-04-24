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

    dp_y <- afscOM::subset_dem_params(om$dem_params, y, d=1, drop=FALSE)

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
            weights = weights,
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
om_selex <- array(NA, dim=c(nyears, 30, 2, 4, nsims))
em_selex <- array(NA, dim=c(nyears, 30, 2, 4, nsims))

for(s in 1:nsims){
    set.seed(seeds[s])
    model_runs <- run_mse(
        om=om,
        mp=mp,
        mse_options=mse_options,
        seed=seeds[s]
    )

    rISS <- get_rISS(nyears, nsamples=10, weight_type = 1, om, model_runs)

    model_dimensions <- afscOM::get_model_dimensions(om$dem_params$sel)

    aggregated_survey_obs <- aggregate_indices(model_runs$survey_obs, nyears)
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

    om_selex[,,,1:2,s] <- om$dem_params$sel[1:nyears,,,1,]
    om_selex[,,,3:4,s] <- om$dem_params$surv_sel[1:nyears,,,1,]

    em_selex[,,,1:2,s] <- sabie_rtmb_model$rep$fish_sel[1,,,,]
    em_selex[,,,3:4,s] <- sabie_rtmb_model$rep$srv_sel[1,,,,]

}

sd <- RTMB::sdreport(sabie_rtmb_model)



get_rISS <- function(nyears, nsamples, weight_type, om, model_runs, obs){
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

                ac <- array(age_comp_obs[y,,1:2,,f,s], c(30*2, 1))
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


om_recs <- apply(model_runs$naa[,1,,,drop=FALSE], 1, sum)
plot(om_recs, type="l", ylim=c(0, 200))
lines(sabie_rtmb_model$rep$Rec[1,], col="red")

om_ssb <- apply(model_runs$naa[1:nyears,,1,,drop=FALSE]*om$dem_params$mat[1:nyears,,1,,drop=FALSE]*om$dem_params$waa[1:nyears,,1,,drop=FALSE], 1, sum)
plot(om_ssb,type="l", ylim=c(0, 300))
lines(sabie_rtmb_model$rep$SSB[1,], col="red")

format_selex <- function(s){
    reshape2::melt(s) %>% as_tibble() %>%
        rename(
            "time"="Var1",
            "age"="Var2",
            "sex"="Var3",
            "fleet"="Var4",
            "sim"="Var5",
            "selex"="value"
        ) %>%
        mutate(
            sex=ifelse(sex==1,"Female","Male"),
            fleet=case_when(
                fleet == 1 ~ "Fixed Fishery",
                fleet == 2 ~ "Trawl Fishery",
                fleet == 3 ~ "Longline Survey",
                fleet == 4 ~ "Trawl Survey"
            )
        ) %>%
        filter(time %in% c(1, 65)) %>%
        mutate(
            timeblock = factor(time, labels=c("Early", "Late"))
        )
}

format_selex(om_selex) %>%
    left_join(format_selex(em_selex), by=c("time", "age", "sex", "fleet", "sim", "timeblock"), suffix=c(".om", ".em")) %>%
    group_by(timeblock, age, sex, fleet) %>%
    median_qi(selex.om, selex.em, .width=c(0.50, 0.85)) %>%

    ggplot()+
        geom_line(aes(x=age, y=selex.om, linetype=timeblock), color="red")+
        geom_lineribbon(aes(x=age, y=selex.em, ymin=selex.em.lower, ymax=selex.em.upper, group=timeblock), size=0.1)+
        scale_fill_brewer(palette="Blues")+
        facet_grid(rows=vars(fleet), cols=vars(sex))+
        custom_theme


est_selex_df <- reshape2::melt(sabie_rtmb_model$rep$srv_sel) %>% as_tibble() %>%
    rename(
        "region"="Var1",
        "time"="Var2",
        "age"="Var3",
        "sex"="Var4",
        "fleet"="Var5"
    ) %>%
    mutate(
        age=age+1,
        sex = factor(ifelse(sex==1, "F", "M")),
        fleet = factor(ifelse(fleet==1, "Fixed", "Trawl")),
        region=as.double(region)
    )


om_selex_df <- reshape2::melt(om$dem_params$surv_sel[1:65,,,1,,drop=FALSE]) %>% as_tibble() %>% mutate(region=1)


selex_df <- est_selex_df %>% left_join(om_selex_df, by=c("time", "age", "sex", "region", "fleet"), suffix=c(".em", ".om"))

ggplot(selex_df %>% filter(time %in% c(1, 65)))+
    geom_line(aes(x=age, y=value.em, group=interaction(sex, time), linetype=factor(time)))+
    geom_line(aes(x=age, y=value.om, group=interaction(sex, time), linetype=factor(time)), color="red")+
    facet_grid(~fleet+sex)+
    custom_theme








lapply(sabie_rtmb_model$rep[grepl("nLL", names(sabie_rtmb_model$rep))], \(x) sum(x))

om_catch <- apply(model_runs$caa, c(1, 5), sum)
PredCatch <- sabie_rtmb_model$rep$PredCatch[1,,]
dimnames(PredCatch) = list("time"=1:82, "fleet"=c("Fixed", "Trawl"))

melt(om_catch) %>% as_tibble() %>% left_join(melt(PredCatch), by=c("time", "fleet")) %>%
    ggplot(aes(x=time))+
        geom_line(aes(y=value.x), color="red")+
        geom_line(aes(y=value.y))+
        facet_wrap(~fleet)+
        coord_cartesian(ylim=c(0, 75))+
        custom_theme


sabie_rtmb_model$rep$NAA[,,,2]

plot(1:30, om$dem_params$waa[1,,1,1], type="l")
lines(1:30, om$dem_params$waa[1,,2,1], col="red")


spock_rps <- SPoCK::Get_Reference_Points(
    data,
    sabie_rtmb_model$rep,
    SPR_x=0.40
)

dp_y <- subset_dem_params(om$dem_params, 65, 1, drop=FALSE)

sel_est <- array(NA, dim=c(82, 30, 2, 1, 2))
sel_est[,,,,1] <- aperm(sabie_rtmb_model$rep$fish_sel, c(2, 3, 4, 1, 5))[,,,,1]
sel_est[,,,,2] <- aperm(sabie_rtmb_model$rep$fish_sel, c(2, 3, 4, 1, 5))[,,,,2]

calculate_joint_selret(
    subset_matrix(sel_est, 65, 1, drop=FALSE),
    subset_matrix(om$dem_params$ret[,,,1,,drop=FALSE], 65, 1, drop=FALSE),
    prop_fs <- c(0.80, 0.20)
)

om_rps <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,1],
    mat = dp_y$mat[,,1,1],
    waa = dp_y$waa[,,1,1],
    sel =  joint_selret$sel[,,1,1,drop=FALSE],
    ret = joint_selret$ret[,,1,1,drop=FALSE],
    avg_rec = mean(apply(model_runs$naa[,1,,], 1, sum))/2,
    spr_target = 0.40
)

om_rps
spock_rps$rep$F_x
spock_rps$rep$SB_F_x*mean(sabie_rtmb_model$rep$Rec)/2
spock_rps$rep$SB0*mean(sabie_rtmb_model$rep$Rec)/2

##### Males NAA
male_NAA_em <- sabie_rtmb_model$rep$NAA[1,,,2]
female_NAA_em <- sabie_rtmb_model$rep$NAA[1,,,1]
male_NAA_om <- apply(model_runs$naa[,,2,], c(1, 2), sum)
female_NAA_om <- apply(model_runs$naa[,,1,], c(1, 2), sum)
dimnames(male_NAA_em) <- list("time"=1:83, "age"=2:31)

reshape2::melt(male_NAA_em) %>%
    left_join(melt(male_NAA_om), by=c("time", "age")) %>%
    ggplot(aes(x=time))+
        geom_line(aes(y=value.x))+
        geom_line(aes(y=value.y), color="red")+
        facet_wrap(~age, scales="free_y")+
        custom_theme

srvtwl_selex_males_om <- male_NAA_om*subset_matrix(om$dem_params$surv_sel, 1, d=4, drop=TRUE)[1:83,,2,2]
srvtwl_selex_males_em_true <- male_NAA_em*subset_matrix(om$dem_params$surv_sel, 1, d=4, drop=TRUE)[1:83,,2,2]
srvtwl_selex_males_em <- male_NAA_em[1:82,]*sabie_rtmb_model$rep$srv_sel[1,,,2,2]#*subset_matrix(om$dem_params$surv_sel, 1, d=4, drop=TRUE)[1:83,,2,2] # sabie_rtmb_model$rep$srv_sel[1,,,2,2]
dimnames(srvtwl_selex_males_em) <- list("time"=1:82, "age"=2:31)
dimnames(srvtwl_selex_males_em_true) <- list("time"=1:83, "age"=2:31)

reshape2::melt(srvtwl_selex_males_em) %>%
    left_join(melt(srvtwl_selex_males_om), by=c("time", "age")) %>%
    left_join(melt(srvtwl_selex_males_em_true), by=c("time", "age")) %>%
    ggplot(aes(x=time))+
        geom_line(aes(y=value.x))+
        geom_line(aes(y=value.y), color="red")+
        geom_line(aes(y=value), color="blue")+
        facet_wrap(~age, scales="free_y")+
        custom_theme


srvtwl_selex_males_em/male_NAA_em







###### OSA RESIDUALS
names(sabie_rtmb_model$rep)
## Not run: 
comp_props <- get_comp_prop(
    data = data, 
    rep = sabie_rtmb_model$rep, 
    age_labels = 2:31, 
    len_labels = seq(41, 99, 2), 
    year_labels = 1960:(1960+nyears-1)
)

osa_results <- get_osa(obs_mat = comp_props$Obs_SrvAge_mat,
                      exp_mat = comp_props$Pred_SrvAge_mat,
                      N = 50,
                      years = 1960:(1960+nyears-1),
                      fleet = 2,
                      bins = 2:31,
                      comp_type = 2,
                      bin_label = "Age")

plot_resids(osa_results)


apply(data$ObsSrvAgeComps[1,,,2,2], 2, mean) # Males
apply(data$ObsSrvAgeComps[1,,,1,2], 2, mean) # Females
