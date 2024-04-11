
rm(list=ls())

library(TMB)
library(tidyverse)
library(ggdist)
setwd("~/Desktop/Projects/SablefishMSE")
source("R/run_mse_multiple.R")
source("R/run_mse.R")
source("R/format_em_data.R")
source("R/reference_points.R")
source("R/harvest_control_rules.R")
source("R/simulate_TAC.R")
source("R/data_utils.R")

library(devtools)
devtools::load_all("~/Desktop/Projects/afscOM")

sable_om <- readRDS("data/sablefish_om.RDS") # Read this saved OM from a file
sable_om$model_options$simulate_observations <- TRUE # Enable simulating observations
# Generate age comp samples as integers rather than proportions
# (this is necesarry for the SpatialSablefishAssessment TMB model)
sable_om$model_options$obs_pars$surv_ll$as_integers = F
sable_om$model_options$obs_pars$surv_tw$as_integers = F
sable_om$model_options$obs_pars$fish_fx$as_integers = F
# Turn off estimation model
sable_om$model_options$run_estimation = FALSE

tier3 <- function(ref_pts, naa, dem_params){
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
    )
}

nsims <- 20
set.seed(1120)
seeds <- sample(1:1e3, nsims)
nyears <- 130
sable_om$model_options$obs_pars$fish_fx$ac_samps <- 30
sable_om$model_options$obs_pars$surv_ll$ac_samps <- 30
sable_om$model_options$obs_pars$surv_tw$ac_samps <- 30
sable_om$model_options$obs_pars$surv_ll$rpn_cv <- 0.10
sable_om$model_options$obs_pars$surv_tw$rpw_cv <- 0.10

s <- 1
mse_tier3 <- run_mse_multiple(nsims=nsims, seeds=seeds, nyears=nyears, om=sable_om, hcr=tier3)

data <- data.frame()

for(s in 1:nsims){
    print(s)
    #' Generate assessment inputs for simulation s
    i=nyears-63
    y <- 63+i
    assess_inputs <- format_em_data(
                            nyears = y,
                            dem_params = afscOM::subset_dem_params(sable_om$dem_params, 64:y, d=1, drop=FALSE),
                            land_caa = afscOM::subset_matrix(mse_tier3$land_caa[64:y,,,,,s,drop=FALSE], 1, d=6, drop=TRUE),
                            survey_indices = afscOM::subset_dem_params(afscOM::subset_dem_params(mse_tier3$survey_obs, 64:y, d=1, drop=FALSE), s, d=5, drop=TRUE),
                            fxfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$fxfish_acs[63:(y-1),,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                            ll_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$ll_acs[63:(y-1),,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                            model_options = sable_om$model_options,
                            added_years = i,
                            file_suffix = s
                        )   
    file.remove(paste0("data/sablefish_em_data_curr_",s,".RDS"))
    file.remove(paste0("data/sablefish_em_par_curr_",s,".RDS"))

    #' Add bias adjustment parameters to data inputs
    assess_inputs$new_data$bias_ramp <- 1
    assess_inputs$new_data$sigma_R <- 1.04
    assess_inputs$new_data$sigma_R_early <- 0.40
    assess_inputs$new_data$sigma_R_early_end <- 15
    assess_inputs$new_data$b_year_start_idx <- 20
    assess_inputs$new_data$b_year_end_idx <- 30
    assess_inputs$new_data$end_yr_rec_est_idx <- length(assess_inputs$new_data$years)-1
    assess_inputs$new_data$bmax <- 1
    assess_inputs$new_data$ba50 <- 5

    #' Do intermediate recruitment calcs
    # bias_ramp <- rep(0, nyears)
    # for(y in (1930-30+2):(2022+nyears-64)){
    #     i <- y - (1930-30+2) +1
    #     if(y < 1980){
    #         bias_ramp[i] <- 0
    #     }else if(y >= 1980 & y < 1990){
    #         bias_ramp[i] <- 1*(1-((y-1980)/(1990-1980)))
    #     }else if(y >= 1990 & y < (1960+nyears-5)){
    #         bias_ramp[i] <- 1
    #     }else if(y >= (1960+nyears-5)){
    #         bias_ramp[i] <- 1*(((y-2022+nyears)/(2022+nyears-(1960+nyears-5))))
    #     }
    # }

    # bias_ramp <- bias_ramp[(length(bias_ramp)-nyears+1):length(bias_ramp)]
    # sigrs <- c(rep(0.40, 15), rep(1.20, 5), rep(1.04, nyears-15-5))
    # recruitment <- apply(mse_tier3$naa[1:nyears,1,,,s], 1, sum)
    # ln_mean_rec <- mean(log(recruitment))
    # ln_rec_devs <- log(recruitment) - ln_mean_rec + bias_ramp*sigrs^2/2

    # ln_init_rec_devs <- c(-0.0201727896823, -0.0211653628442, -0.0222762796367, -0.0235167887645, -0.0248376531671, -0.0263140331771, -0.0279775591988, -0.0298335784895, -0.0318628708855, -0.0341045755794, -0.0365927819101, -0.0393071338025, -0.0422739844026, -0.0455070747877, -0.0490012744867, -0.0527273187832, -0.0566564781594, -0.0607215085617, -0.0647839822855, -0.0686725000884, -0.0721473710816, -0.0747356469988, -0.0756711664860, -0.0736932654205, -0.0674020989957, -0.0545598257407, -0.0326432856417, 0.000831967151461, 0.0477478898283, 0.106462820860, 0.170899944240, 0.228214073308, 0.260537688760, 0.272295671697, 0.223888837810, 0.138327569143, 0.0409977654323, -0.0504299349341, -0.139913190910, -0.225833906653, -0.298717396287, -0.351735855526, -0.363169725158, -1.08102049603, -0.951503149105, -0.900354062579, -0.671126924421, 1.14538462882, 0.993817636166, 0.167556465633, 1.45911595657, 0.686962952511, -0.296119130703, -0.240150999009, 0.105223647861, -0.806134517718, -1.20159915722, -1.10370002507, -0.280362194813, 0.385861443905, -0.677219759900, 0.451877334609, -0.816981612496, -0.754631289491, -0.241530322750, 0.269235393826, -0.431497255133, 0.853730546539, 0.0317112429988, 0.0425907663835, 1.01280894948, -0.151876454291, -0.447976422128, -0.272519813754, -0.712160267655, -0.413806228641, -0.468138784792, 0.0159327161347, 0.319484365520, -0.407168000172, -0.321847602343, -1.19277559360, -0.726547291688, -0.0723645963570, 1.16309263388, 0.337204769385, 1.81806716488, 1.85659954644, 1.23682827751, 1.98814135583, 0.468638127271)[1:28]
    # ass_rec <- exp(3.28521376476 + ln_init_rec_devs - 0.40^2/2)
    # scaled_ln_init_rec_devs <- log(ass_rec) - ln_mean_rec + 0.40^2/2

    #' Do intermediate F dev calculations
    # true_fs_ll <- apply(mse_tier3$faa[1:nyears,,,1,1,s], 1, max)
    # ln_ll_F_avg <- log(mean(true_fs_ll))
    # ln_ll_F_devs <- log(true_fs_ll) - ln_ll_F_avg

    # true_fs_tw <- apply(mse_tier3$faa[1:nyears,,,1,2,s], 1, max)
    # ln_tw_F_avg <- log(mean(true_fs_tw))
    # ln_tw_F_devs <- log(true_fs_tw) - ln_tw_F_avg

    #' Map off parameters
    par_map <- list(
        # Turn off japanses fishery/survey and CPUE parameters
        ln_srv_jap_ll_sel_pars = factor(rep(NA, length(assess_inputs$new_parameters$ln_srv_jap_ll_sel_pars))),
        ln_srv_jap_fishery_ll_sel_pars = factor(rep(NA, length(assess_inputs$new_parameters$ln_srv_jap_fishery_ll_sel_pars))),
        ln_srv_jap_ll_q = factor(rep(NA, length(assess_inputs$new_parameters$ln_srv_jap_ll_q))),
        ln_srv_jap_fishery_ll_q = factor(rep(NA, length(assess_inputs$new_parameters$ln_srv_jap_fishery_ll_q))),
        ln_ll_cpue_q = factor(rep(NA, length(assess_inputs$new_parameters$ln_ll_cpue_q))),

        # Turn off recruitment related parameters
        #ln_mean_rec     = factor(NA),
        #ln_ll_F_avg     = factor(NA),
        #ln_rec_dev      = factor(rep(NA, length(ln_rec_devs))),
        # ln_init_rec_dev = factor(rep(NA, length(ln_init_rec_devs))),

        # Turn off F deviates
        # ln_ll_F_avg     = factor(NA),
        # ln_trwl_F_avg   = factor(NA),
        # ln_ll_F_devs    = factor(rep(NA, length(assess_inputs$new_parameters$ln_ll_F_devs))),
        # ln_trwl_F_devs    = factor(rep(NA, length(assess_inputs$new_parameters$ln_trwl_F_devs))),

        # Turn off LL selectivities
        #ln_ll_sel_pars = factor(c(1:9, 4, 10, 11)), # share delta for first time block across sexes
        ln_trwl_sel_pars = factor(c(1,2,3,2)), # share power parameter across sexes
        ln_ll_sel_pars = factor(rep(NA, length(assess_inputs$new_parameters$ln_ll_sel_pars))),
        # Turn off mortality related parameters
        # ln_M            = factor(NA),
        ln_M_year_devs  = factor(rep(NA, length(assess_inputs$new_parameters$ln_M_year_devs))),
        ln_M_age_devs   = factor(rep(NA, length(assess_inputs$new_parameters$ln_M_age_devs)))
    )
    
    # fishery ll selx
    # sel_pars <- array(c(1.4257e+000, 1.2222e+000, 6.5829e-001,
    #                     -6.9100e-001, 5.3936e-001 , 8.1499e-001,
    #                     1.9380e+000, 1.4919e+000, 1.0631e+000,
    #                     -6.9100e-001, -9.7358e-002, -3.7957e-001), dim = c(3,2,2))

    #assess_inputs$new_parameters$ln_mean_rec <- ln_mean_rec
    #assess_inputs$new_parameters$ln_ll_F_avg <- -3.0236
    #assess_inputs$new_parameters$ln_init_rec_dev <- scaled_ln_init_rec_devs
    #assess_inputs$new_parameters$ln_M <- log(sable_om$dem_params$mort[1,1,1,1])
    #assess_inputs$new_parameters$ln_rec_dev <- ln_rec_devs
    # assess_inputs$new_parameters$ln_ll_sel_pars <- sel_pars
    # assess_inputs$new_parameters$ln_ll_F_avg <- ln_ll_F_avg
    # assess_inputs$new_parameters$ln_trwl_F_avg <- ln_tw_F_avg
    # assess_inputs$new_parameters$ln_ll_F_devs <- ln_ll_F_devs
    # assess_inputs$new_parameters$ln_trwl_F_devs <- ln_tw_F_devs
    # assess_inputs$new_data$loglik_wgt_ll_catch <- 100
    # assess_inputs$new_data$loglik_wgt_trwl_catch <- 100
    assess_inputs$new_data$trwl_sel_type <- 1
    assess_inputs$new_parameters$ln_trwl_sel_pars <- array(c(2.10484881040, 0.00529011908453,
                                                             1.77013905300, 0.00529011908453), dim=c(1, 2, 2))

    assess_inputs$new_data$loglik_wgt_ll_catchatage <- 0.797868466479416
    assess_inputs$new_data$loglik_wgt_srv_dom_ll_age <- 3.72382880611179
    assess_inputs$new_data$loglik_wgt_srv_jap_ll_age <- 1.27151115308662
    assess_inputs$new_data$loglik_wgt_ll_catchatlgth_m <- 5.21629903889062
    assess_inputs$new_data$loglik_wgt_ll_catchatlgth_f <- 4.94489032547033
    assess_inputs$new_data$loglik_wgt_srv_dom_ll_lgth_m <- 1.11519421735366
    assess_inputs$new_data$loglik_wgt_srv_dom_ll_lgth_f <- 1.49981566417443
    assess_inputs$new_data$loglik_wgt_srv_jap_ll_lgth_m <- 0.902165573415091
    assess_inputs$new_data$loglik_wgt_srv_jap_ll_lgth_f <- 1.26841555181968
    assess_inputs$new_data$loglik_wgt_trwl_catchatlgth_m <- 0.255150139457704
    assess_inputs$new_data$loglik_wgt_trwl_catchatlgth_f <- 0.35007890900583
    assess_inputs$new_data$loglik_wgt_srv_nmfs_trwl_lgth_m <- 0.450498802823617
    assess_inputs$new_data$loglik_wgt_srv_nmfs_trwl_lgth_f <- 0.672532119442201

    #' Fit the model
    compile("inst/CurrentAssessment.cpp")
    dyn.load(dynlib("inst/CurrentAssessment"))

    my_model = TMB::MakeADFun(data = assess_inputs$new_data,
                              parameters = assess_inputs$new_parameters,
                              map = par_map,
                              DLL = "CurrentAssessment",
                              silent = TRUE)
    mle_optim = tryCatch(
         expr = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000)),
         error = function(e){e}
    )

    # if(inherits(mle_optim, "error")){
    #     next;
    # }
    
    # tryCatch(expr =
    #     for(i in 1:5) {
    #         g = as.numeric(my_model$gr(mle_optim$par))
    #         h = optimHess(mle_optim$par, fn = my_model$fn, gr = my_model$gr)
    #         mle_optim$par = mle_optim$par - solve(h,g)
    #         mle_optim$objective = my_model$fn(mle_optim$par)
    #     }
    # , error = function(e){e})

    sdr <- sdreport(my_model)
    if(any(is.nan(sdr$sd))){
        print("Model convergence failed")
        next;
    }
    
    mle_report = my_model$report(my_model$env$last.par.best)
    
    # plot(mle_report$sel_ll_f[,3], ylim = c(0,1), type = 'l')
    # lines(sable_om$dem_params$sel[164,,1,1,1], ylim = c(0,1))
    
    # plot(mle_report$sel_ll_m[,3], ylim = c(0,1), type = "l")
    # lines(sable_om$dem_params$sel[164,,2,1,1], ylim = c(0,1))
    
    # plot(mle_report$sel_trwl_f[,1], ylim = c(0,1))
    # lines(sable_om$dem_params$sel[1,,1,1,2], ylim = c(0,1))
    
    # plot(mle_report$sel_trwl_m[,1], ylim = c(0,1))
    # lines(sable_om$dem_params$sel[1,,2,1,2], ylim = c(0,1))

    est_rec <- SpatialSablefishAssessment::get_recruitment(mle_report)
    om_rec <- apply(mse_tier3$naa[1:nyears,1,,,s], 1, sum)

    est_ssb <- SpatialSablefishAssessment::get_SSB(mle_report)
    om_ssb <- apply(mse_tier3$naa[1:nyears,,1,1,s]*sable_om$dem_params$mat[1:nyears,,1,1]*sable_om$dem_params$waa[1:nyears,,1,1], 1, sum)

    est_Fs <- SpatialSablefishAssessment::get_fishing_mortalities(mle_report) %>% as_tibble() %>% group_by(Year) %>% summarise(F = sum(F))
    om_Fs <- mse_tier3$out_f[1:nyears,1, 1, 1, s]

    est_fish_Fs <- SpatialSablefishAssessment::get_fishing_mortalities(mle_report) %>% as_tibble() %>%
        mutate(
            true_Fs = c(apply(mse_tier3$faa[1:nyears,,,1,1,s], 1, max), apply(mse_tier3$faa[1:nyears,,,1,2,s], 1, max)),
            avg_F = c(rep(exp(mle_optim$par["ln_ll_F_avg"]), nyears), rep(exp(mle_optim$par["ln_trwl_F_avg"]), nyears)),
            true_avg_F = c(rep(exp(ln_ll_F_avg), nyears), rep(exp(ln_tw_F_avg), nyears))
        )

    est_catches <- SpatialSablefishAssessment::get_catches(mle_report) %>% as_tibble() %>%
        pivot_wider(names_from=type, values_from=Catch)


    tmp <- data.frame(
        Year=1960:(1960+nyears-1), 
        est_rec=est_rec$Recruitment, 
        om_rec=om_rec, 
        est_ssb=est_ssb$SSB, 
        om_ssb=om_ssb, 
        est_F=est_Fs$F, 
        om_F=om_Fs, 
        est_ll_f=est_fish_Fs %>% filter(Fishery == "Fixed gear") %>% pull(F),
        om_ll_f=est_fish_Fs %>% filter(Fishery == "Fixed gear") %>% pull(true_Fs),
        est_tw_f=est_fish_Fs %>% filter(Fishery == "Trawl gear") %>% pull(F),
        om_tw_f=est_fish_Fs %>% filter(Fishery == "Trawl gear") %>% pull(true_Fs),
        est_ll_catch = est_catches %>% filter(Fishery == "Fixed gear") %>% pull(Predicted),
        om_ll_catch = est_catches %>% filter(Fishery == "Fixed gear") %>% pull(Observed),
        est_tw_catch = est_catches %>% filter(Fishery == "Trawl") %>% pull(Predicted),
        om_tw_catch = est_catches %>% filter(Fishery == "Trawl") %>% pull(Observed),
        sim=s)

    data <- bind_rows(data, tmp)

}


df <- data %>% group_by(Year) %>%
    ggdist::median_qi(est_rec, om_rec, est_ssb, om_ssb, est_F, om_F, est_ll_f, om_ll_f, est_tw_f, om_tw_f, est_ll_catch, om_ll_catch, est_tw_catch, om_tw_catch, .width=c(0.50, 0.80))

rec_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_rec, ymin=om_rec.lower, ymax=om_rec.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_rec, ymin=est_rec.lower, ymax=est_rec.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 200))+
    ggtitle("Recruitment")+
    guides(fill="none")+
    theme_bw()

ssb_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_ssb, ymin=om_ssb.lower, ymax=om_ssb.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_ssb, ymin=est_ssb.lower, ymax=est_ssb.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 350))+
    ggtitle("Spawning Biomass")+
    guides(fill="none")+
    theme_bw()

F_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_F, ymin=om_F.lower, ymax=om_F.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_F, ymin=est_F.lower, ymax=est_F.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 0.3))+
    ggtitle("Fishing Mortality")+
    guides(fill="none")+
    theme_bw()

ll_F_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_ll_f, ymin=om_ll_f.lower, ymax=om_ll_f.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_ll_f, ymin=est_ll_f.lower, ymax=est_ll_f.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 0.12))+
    ggtitle("Fixed Gear Fishing Mortality")+
    guides(fill="none")+
    theme_bw()

tw_F_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_tw_f, ymin=om_tw_f.lower, ymax=om_tw_f.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_tw_f, ymin=est_tw_f.lower, ymax=est_tw_f.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 0.12))+
    ggtitle("Trawl Gear Fishing Mortality")+
    guides(fill="none")+
    theme_bw()

ll_catch_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_ll_catch, ymin=om_ll_catch.lower, ymax=om_ll_catch.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_ll_catch, ymin=est_ll_catch.lower, ymax=est_ll_catch.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 50))+
    ggtitle("Fixed Gear Catch")+
    guides(fill="none")+
    theme_bw()

tw_catch_plot <- ggplot(df)+
    geom_lineribbon(aes(x=Year, y=om_tw_catch, ymin=om_tw_catch.lower, ymax=om_tw_catch.upper), size=0.6)+
    geom_line(aes(x=Year, y=est_tw_catch, ymin=est_tw_catch.lower, ymax=est_tw_catch.upper), color="red")+
    scale_fill_brewer()+
    scale_y_continuous(limits=c(0, 50))+
    ggtitle("Trawl Gear Catch")+
    guides(fill="none")+
    theme_bw()

library(patchwork)

library(cowplot)

plot_grid(
    plot_grid(rec_plot, ssb_plot, F_plot, ncol=3),
    plot_grid(
        plot_grid(ll_F_plot, tw_F_plot, ncol=1), 
        plot_grid(ll_catch_plot, tw_catch_plot, ncol=1), 
        ncol=2
    ),
    ncol=1
)

ggsave("~/Desktop/em_bias_summary1.pdf", width=11, unit="in")


# par(mfrow=c(1, 2))
# model_sels <- my_model$report()$sel_ll_f
# om_sels <- t(sable_om$dem_params$sel[c(1, 40, 60), , 1, 1, 1])
# plot(1:30, model_sels[,3], linetype="dashed", type="l", ylim=c(0, 1), main="Female LL Selectivity")
# lines(1:30, om_sels[,3], col="red", linetype="dashed")
# plot(1:30, my_model$report()$sel_ll_m[,3], linetype="dashed", type="l", ylim=c(0, 1), main="Male LL Selectivity")
# lines(1:30, t(sable_om$dem_params$sel[c(1, 40, 60), , 2, 1, 1])[,3], col="red", linetype="dashed")

# model_sels <- my_model$report()$sel_trwl_f
# om_sels <- t(sable_om$dem_params$sel[c(1), , 1, 1, 2])
# plot(1:30, model_sels[,1], type="l", ylim=c(0, 1))
# lines(1:30, om_sels, col="red", linetype="dashed")

# t(sable_om$dem_params$sel[c(1), , 1, 1, 2])
# sel_50 <- exp(2.013419)
# delta <- exp(2.215031)
# ages <- 1:31
# ((ages/sel_50)^(sel_50/(0.5*sqrt(sel_50^2 + 4*delta^2) - sel_50)))*exp((sel_50-ages)/(0.5*(sqrt(sel_50^2+4*delta^2)-sel_50)))

# (pow(ages[age] / sel_50, sel_50 / (0.5*(sqrt(square(sel_50)+4*square(delta))-sel_50)))*exp((sel_50-ages[age])/(0.5*(sqrt(square(sel_50)+4*square(delta))-sel_50))));

# names(my_model$report())[grepl("catch", names(my_model$report()))]

# my_model$report()$annual_ll_catch_pred
# my_model$report()$ll_fishery_catch
# plot(1:nyears, my_model$report()$ll_fishery_catch, type="l")
# lines(1:nyears, my_model$report()$annual_ll_catch_pred, col="red")

# my_model$report()$annual_trwl_catch_pred
# my_model$report()$trwl_fishery_catch
# plot(1:nyears, my_model$report()$trwl_fishery_catch, type="l")
# lines(1:nyears, my_model$report()$annual_trwl_catch_pred, col="red")

# 1/(1+exp(-(exp(0.8149))*(1:31-exp(0.65829))))

# sdreport(my_model)

# my_model$par[grepl("ln_ll_sel_pars", names(my_model$par))]
# assess_inputs$new_parameters$ln_ll_sel_pars

# my_model$report()$sel_trwl_f[,1]
# # Females block 1
# 1/(1+exp(-(exp(-0.690999985510))*(1:31-exp(1.42565337612))))
# # Femaes block 2
# 1/(1+exp(-(exp(0.539364859399))*(1:31-exp(1.22222290063))))
# # Females block 3
# 1/(1+exp(-(exp(0.8149))*(1:31-exp(0.65829))))

# sel_pars <- array(NA, dim=c(3, 2, 2))
# sel_pars[1,1,2] <- 1.42565337612
# sel_pars[1,2,2] <- -0.690999985510
# sel_pars[2,1,2] <- 1.22222290063
# sel_pars[2,2,2] <- 0.539364859399
# sel_pars[3,1,2] <- 0.65829
# sel_pars[3,2,2] <- 0.8149

# sel_pars[1,1,1] <- 1.93799670819
# sel_pars[1,2,1] <- 0.000808317638002
# sel_pars[2,1,1] <- 1.49186769564
# sel_pars[2,2,1] <- -0.0973584565150
# sel_pars[3,1,1] <- 1.06312413114
# sel_pars[3,2,1] <- -0.379568217967

