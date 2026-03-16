fs <- list.files(file.path(here::here(), "data", "active"), full.names = TRUE)
if(!is.null(hcr_filter))
    fs <- unlist(sapply(hcr_filter, \(x) fs[grepl(paste0(sub("|", "\\|", sub("/", "", sub(" ", "_", tolower(x))), fixed=TRUE), "_\\d+"), fs)]))

o <- tibble::tibble()

k <- 1
for(i in seq_along(fs)){
    x <- fs[i]
    print(x)
    m <- readRDS(x)
    mse <- m$mse_objects
    model_run <- list(mse[[length(mse)]])[[1]]
    extra_columns <- expand.grid(om=model_run$om$name, hcr=model_run$mp$name)

    seeds <- m$seeds

    convergence <- model_run$converged

    fullconvergence <- reshape2::melt(convergence, value.name="converged") %>% as_tibble() %>%
        rename(sim_year="Var1", sim="Var2") %>%
        mutate(om=extra_columns$om, hcr=extra_columns$hcr) %>%
        mutate(sim = factor(sim, labels=seeds))
        # group_by(om, hcr, sim) %>%
        # summarise(n_converged=sum(converged))
    
    o <- bind_rows(o, fullconvergence)

    if(i %% 500 == 0){
        o %>% write_csv(file=file.path(here::here(), "data", paste0("zahneretal_26_spatialmse_fluxes_", k, ".csv")))
        k <- k + 1
        o <- tibble::tibble()
    }

}

convergence <- read_csv(file.path("data", "zahneretal_26_spatialmse_convergence.csv"))

total_sims <- 51*200*8*6 # 51 years, 200 sims, 8 hcrs, 6 oms

convergence %>% filter(!converged) %>% nrow() / total_sims

unconverged <- convergence %>% filter(!converged) # 0.01527% failure rate

unconverged %>% group_by(om, hcr) %>% 
    summarise(n_unconverged=n()) %>% 
    arrange(desc(n_unconverged)) %>%
    print(n=50)

unconverged %>% distinct(sim) # 67 sims have at least one convergence failure; 33.5%

unconverged %>% group_by(sim, om, hcr) %>%
    summarise(n_unconverged=n()) %>%
    arrange(desc(n_unconverged)) %>%
    print(n=50)

# seeds with >= 10 unconverged EM runs: 76, 1041, 1583, 1027, 1796
# seeds with >= 1% unconverged EM runs: 76, 1041, 1583, 1027
unconverged %>% group_by(sim) %>%
    summarise(n_unconverged=n()) %>%
    mutate(perc = n_unconverged / (51*8*6)) %>%
    arrange(desc(perc)) %>%
    print(n=50)

# First 10 years most likely to fail convergence diagnostics
unconverged %>% count(sim_year) %>% arrange(sim_year) %>% print(n=51)

unconverged_diags <- get_unconverged_model_diags(unconverged, om_list, seed_list)

# 1 "Unknown" non converged = Non finite SE (sim 627)
# 596 non converged bc of parameter correlations
# 151 non converged bc of absolute gradient size 
unconverged_diags %>% count(parameter)

unconverged_diags %>% filter(
    parameter %in% 
    c("ln_F_mean", "ln_M", "ln_global_R0", "ln_fish_fixed_sel_pars", "ln_srv_fixed_sel_pars")
) %>% print(n=200)

unconverged_diags %>% filter(
    parameter %in%
    c("ln_srv_fixed_sel_pars and ln_srv_fixed_sel_pars")
)

# Four models have absolute gradient components >0.01 (811, 580, 757, 1108)
# Simulations 811 ad 580 have gradient components >1 and should be removed from analysis
unconverged_diags %>% filter(
        parameter %in% 
        c("ln_F_mean", "ln_M", "ln_global_R0", "ln_fish_fixed_sel_pars", "ln_srv_fixed_sel_pars")
    ) %>% 
    separate(problem, into=c("problem", "value"), sep=" = ") %>%
    arrange(desc(value))

# Generally, all unconverged models have correlations > 0.999
unconverged_diags %>% filter(
    parameter %in%
    c("ln_srv_fixed_sel_pars and ln_srv_fixed_sel_pars")
) %>% 
    separate(problem, into=c("problem", "value"), sep=" of ")

#########################
# Models that need diagnostic plots
#########################
bad_models <- unconverged_diags %>%
    filter(
        parameter %in%
        c("Unknown", "ln_F_mean", "ln_fish_fixed_sel_pars")
    ) %>%
    bind_rows(
        unconverged_diags %>% filter(
            parameter %in%
            c("ln_srv_fixed_sel_pars and ln_srv_fixed_sel_pars")
        ) %>% slice(1)
    ) %>%
    bind_rows(
        unconverged_diags %>% filter(
            parameter %in% 
            c("ln_M", "ln_global_R0", "ln_srv_fixed_sel_pars")
        ) %>% 
        separate(problem, into=c("problem", "value"), sep=" = ") %>%
        arrange(desc(value)) %>%
        slice(c(1, 2, 24))
    )

for(i in 1:nrow(bad_models)){
    p <- plot_diags(
        hcr = as.character(bad_models[i, "hcr"]),
        om = as.character(bad_models[i, "om"]),
        simulation_seed = as.numeric(bad_models[i, "sim"]),
        simulation_year = as.numeric(bad_models[i, "sim_year"]),
        om_list,
        seed_list,
        report_convergence = TRUE
    )
    show(p)
    ggsave(file.path("figures", "diags", paste0("diags_",i,".jpeg")), width=12, height=12, units="in")
}


