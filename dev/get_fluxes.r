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

    out_dims <- dim(model_run$naa)

    nyears <- out_dims[1]
    nsims <- out_dims[5]

    naa_postmove <- array(0, dim=c(5, 5, 2, 30, nyears, nsims))
    naa_postmort <- array(0, dim=c(5, 5, 1, 30, 2))

    for(sim in 1:nsims){
        # print(sim)
        naa <- model_run$naa[,,,,sim]
        faa <- model_run$faa[,,,,,sim]

        for(year in 2:nyears){
            dp <- afscOM::subset_dem_params(
                model_run$om$dem_params, 
                year, d=1, drop=FALSE
            )
            for(r in 1:5){
                dp.r <- afscOM::subset_dem_params(
                    dp, 
                    r, d=4, drop=FALSE
                )

                naa_tmp <- afscOM::simulate_population(
                    prev_naa = naa[year-1,,,r,drop=FALSE],
                    faa = faa[year-1,,,r,,drop=FALSE],
                    recruitment = naa[year,1,,r,drop=FALSE],
                    dem_params = dp.r,
                    options = model_run$om$model_options
                )

                naa_postmort[r,r,,,] <- naa_tmp$naa
            }
            naa_postmove_tmp <- do_movement(naa_postmort, dp$movement)
            naa_postmove[,,,,year-1, sim] <- naa_postmove_tmp
        }
    }


    dimnames(naa_postmove) <- list(
        "region_start"=c("BS", "AI", "WGOA", "CGOA", "EGOA"),
        "region_end"=c("BS", "AI", "WGOA", "CGOA", "EGOA"),
        "sex"=c("F", "M"),
        "age"=2:31,
        "time"=1:nyears,
        "sim"=1:nsims
    )

    out <- reshape2::melt(naa_postmove, value.name="abund_flux") %>% 
                as_tibble() %>% 
                filter_times(time_horizon) %>%
                mutate(
                    class = case_when(
                        age < 7 ~ "Juvenile",
                        age >= 7 & age < 16 ~ "Adult",
                        age >= 16 ~ "Old"
                    ),
                    class = factor(class, levels=c("Juvenile", "Adult", "Old"))
                ) %>%
                left_join(
                    reshape2::melt(model_run$om$dem_params$waa, value.name="waa") %>% as_tibble(),
                    by=c("time", "age", "sex", "region_start"="region"),
                ) %>%
                left_join(
                    reshape2::melt(model_run$om$dem_params$mat, value.name="mat") %>% as_tibble(),
                    by=c("time", "age", "sex", "region_start"="region"),
                ) %>%
                mutate(
                    bio_flux = abund_flux * waa * mat,
                ) %>%
                group_by(time, sim, class, region_start, region_end) %>%
                summarise(
                    abund_flux=sum(abund_flux),
                    bio_flux=sum(bio_flux)
                ) %>%
                group_by(time, sim, class, region_start) %>%
                mutate(
                    abund_flux_in = sum(abund_flux[region_start != region_end]),
                    bio_flux_in = sum(bio_flux[region_start != region_end])
                ) %>%
                group_by(time, sim, class, region_end) %>%
                mutate(
                    abund_flux_out = -sum(abund_flux[region_start != region_end]),
                    bio_flux_out = -sum(bio_flux[region_start != region_end])
                ) %>%
                mutate(
                    net_abund_flux = abund_flux_in + abund_flux_out,
                    net_bio_flux = bio_flux_in + bio_flux_out
                ) %>%
                bind_cols(extra_columns) %>%
                separate(om, c("Recruitment", "Movement"), sep=c("\\s[|]\\s"), remove=FALSE) %>%
                mutate(
                    Recruitment = factor(Recruitment, levels=c("BH Recruit", "Regime Recruit", "Crash Recruit")),
                    Movement = factor(Movement, levels=c("AB Move", "Climate Move"))
                )
    
    o <- bind_rows(o, out)

    # if(i %% 100 == 0){
    #     o %>% write_csv(file=file.path(here::here(), "data", paste0("zahneretal_26_spatialmse_fluxes_", k, ".csv")))
    #     k <- k + 1
    #     o <- tibble::tibble()
    # }

}



fluxes <- read_csv(file.path("data", "zahneretal_2026_spatialmse_fluxesf00.csv"))

out <- fluxes %>%
    mutate(
        class = case_when(
            age < 7 ~ "Juvenile",
            age >= 7 & age < 16 ~ "Adult",
            age >= 16 ~ "Old"
        ),
        class = factor(class, levels=c("Juvenile", "Adult", "Old"))
    ) %>%
    filter(region_start != region_end) %>%
    group_by(region_start, region_end, time, sim, class) %>%
    summarise(flux=sum(flux)) %>%
    group_by(region_start, region_end, time, class) %>%
    median_qi(flux, .width=c(0.50, 0.80))


ggplot(o, aes(x=time, y=bio_flux, fill=class))+
    geom_col()+
    scale_y_continuous("Abundance Exchange (millions)")+
    facet_grid(rows=vars(region_start), cols=vars(region_end))+
    custom_theme
ggsave("~/Desktop/abundance_exchange.jpeg")

f <- fluxes %>% 
    # filter_hcr_om(hcrs=c("F40"), oms=publication_oms) %>%
    # filter(region_start != region_end) %>%
    group_by(time, sim, class, region_start, region_end, Recruitment, Movement, hcr) %>% 
    summarise(flux=sum(flux)) %>%
    group_by(time, sim, class, region_start, Recruitment, Movement, hcr) %>%
    mutate(flux_in = sum(flux)) %>%
    group_by(time, sim, class, region_end, Recruitment, Movement, hcr) %>%
    mutate(flux_out = -sum(flux)) %>%
    group_by(time, class, region_start, Recruitment, Movement, hcr) %>%
    median_qi(flux_in, flux_out, .width=0.50)

exchange <- fluxes %>% select(!contains("abund")) %>%
    filter(region_start != region_end) %>%
    group_by(time, class, region_start, Recruitment, Movement, hcr) %>%
    median_qi(bio_flux_in, bio_flux_out, net_bio_flux, .width=c(0.50))

retained <- fluxes %>% select(!contains("abund")) %>%
    filter(region_start == region_end) %>%
    mutate(retained=bio_flux-net_bio_flux) %>%
    group_by(time, class, region_start, Recruitment, Movement, hcr) %>%
    median_qi(retained, .width=c(0.50))

net <- retained+exchange

df <- data.frame(x=90, y=c(30, -18), label=c("Flux In", "Flux Out"))
ggplot(fluxes)+
    geom_col(aes(x=time, y=bio_flux_in, fill=class))+
    geom_col(aes(x=time, y=bio_flux_out, fill=class), alpha=0.6)+
    # geom_text(aes(x=90, y=110, label="Flux In"))+
    # geom_text(data=df, aes(x=x, y=y, label=label))+
    scale_y_continuous("Absolute Abundance Exchange (millions)")+
    ggh4x::facet_nested(
        rows=vars(region_start),
        cols=vars(Recruitment, Movement)
    )+
    # facet_wrap(~region_start)+
    custom_theme
ggsave("~/Desktop/net_flux.jpeg")

df <- data.frame(x=90, y=c(12, -6.5), label=c("In", "Out"))
ggplot(exchange)+
    geom_col(aes(x=time, y=net_bio_flux, fill=class))+
    # geom_col(data=retained, aes(x=time, y=retained), fill="grey", alpha=0.50)+
    geom_hline(yintercept=0, color="black")+
    # geom_text(data=df, aes(x=x, y=y, label=label), size=6)+
    scale_y_continuous("Net SSB Exchange (1000s mt)")+
    ggh4x::facet_nested(
        rows=vars(region_start),
        cols=vars(Recruitment, Movement)
        # scales="free_y"
    )+
    # ggh4x::facetted_pos_scales(
    #     y = list(
    #         scale_y_continuous(limits=c(-10, 70), breaks=seq(-10, 70, 25)),
    #         scale_y_continuous(limits=c(-10, 70), breaks=seq(-10, 70, 25)),
    #         scale_y_continuous(limits=c(-10, 35), breaks=seq(-10, 35, 15)),
    #         scale_y_continuous(limits=c(-10, 90), breaks=seq(-10, 90, 25)),
    #         scale_y_continuous(limits=c(-10, 70), breaks=seq(-10, 70, 25))
    #     )
    # )+
    custom_theme

ggplot2::last_plot() +
     plot_annotation(title="Net SSB Exchange by Region under F=0 Harvest Policy",
     caption="Net SSB exchange (inflow - outflow) shown in color; retained SSB shown in grey") & 
     theme(legend.position = "bottom")

ggsave(file.path(figures_dir, paste0("net_bio_flux_", filetype)), width=12, height=12, units="in")


exchange %>% group_by(region_start, Recruitment, Movement, hcr) %>%
    median_qi(bio_flux_in, bio_flux_out, .width=c(0.50), .simple_names=FALSE) %>%
    left_join(
        retained %>% group_by(region_start, Recruitment, Movement, hcr) %>%
            median_qi(retained, .width=c(0.50), .simple_names=FALSE),
        by=c("region_start", "Recruitment", "Movement", "hcr", ".width", ".point", ".interval")
    ) %>%
    group_by(region_start, Recruitment, Movement, hcr) %>%
    # summarise(
    #     total = bio_flux_in + retained,
    #     perc_in = bio_flux_in/total,
    #     perc_retained = retained/total,
    # ) %>%
    pivot_longer(c("bio_flux_in", "retained"), names_to="name") %>%
    ggplot()+
        geom_col(aes(x=Movement, y=value, fill=name))+
        scale_y_continuous("SSB (1000s mt)")+
        facet_grid(region_start~Recruitment)+
        custom_theme

ggsave(file.path(figures_dir, paste0("exchange_retained_props", filetype)), width=12, height=12, units="in")




do_movement <- function(naa_postmort, move){
    move[,,1,1,] <- diag(5)
    return(
        vapply(
            1:30, 
            # Apply movement to the ages individualy
            function(a) {
                # naa_postmort[,,1,a] %*% move[,,1,a,1]
                sapply(
                    1:2,
                    # Apply movement to the sexes individually)
                    function(s) naa_postmort[,,1,a,s] %*% move[,,1,a,s]
                )
            } ,
            FUN.VALUE = array(0, dim=c(5, 5, 2))
        )
    )
}

f %>% group_by(time, class, Recruitment, Movement, hcr) %>%
    summarise(v = sum(net_flux))
