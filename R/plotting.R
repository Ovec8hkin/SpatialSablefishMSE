plot_ssb <- function(data){
    # Plot spawning biomass from OM and EM
    d <- data %>%
        # SSB is females only
        filter(sex == "F") %>%
        # summarise SSB across year and sim 
        group_by(time, hcr, sim, L1) %>%
        summarise(spbio=sum(spbio)) %>%
        # Compute quantiles of SSB distribution
        group_by(time, hcr, L1) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=3)


    plot <- ggplot(d %>% filter(L1 == "naa")) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=hcr, color=hcr), size=0.4)+
        geom_pointrange(data = d %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_vline(xintercept=64, linetype="dashed")+
        geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 300))+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)
}

plot_fishing_mortalities <- function(data){
    # Plot fishing mortality rates from OM and EM
    f <- data %>%
        group_by(time, fleet, L1, hcr) %>%
        median_qi(F, total_F, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=4) %>%
        filter(name == "total_F")

    plot <- ggplot(f %>% filter(L1 == "faa")) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=hcr, color=hcr), size=0.4)+
        geom_pointrange(data = f %>% filter(L1 == "faa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_vline(xintercept=64, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)
}

plot_recruitment <- function(data){
    
    r <- data %>%
        # summarise SSB across year and sim 
        group_by(time, hcr, L1) %>%
        median_qi(rec, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=3)

    mean_rec <- r %>% pull(median) %>% mean

    plot <- ggplot(r %>% filter(L1 == "naa")) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=hcr, color=hcr), size=0.4)+
        geom_pointrange(data = r %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_hline(yintercept = mean_rec, linetype="dashed") + 
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 120))+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)
}

plot_phase_diagram <- function(model_runs, extra_columns, dem_params, nyears){
    d <- get_ssb_biomass(model_runs, extra_columns, dem_params = ) %>%
        # SSB is females only
        filter(sex == "F", L1 == "naa") %>%
        # summarise SSB across year and sim 
        group_by(time, hcr, sim, L1) %>%
        summarise(spbio=sum(spbio)) %>%
        left_join(
            reshape2::melt(mse_small$out_f),
            by = c("time", "sim"),
        ) %>%
        group_by(time, hcr) %>%
        median_qi(spbio, value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        drop_na()

    ref_pts <- get_reference_points(model_runs, extra_columns, dem_params, year=nyears)

    plot <- ggplot(d)+
        geom_point(aes(x=spbio, y=value, color=time), size=1.5)+
        geom_segment(
            aes(x=spbio, y=value, xend=after_stat(lead(x)), yend=after_stat(lead(y)), color=time), 
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(yintercept=ref_pts$Fref, linetype="dashed")+
        geom_vline(xintercept=ref_pts$Bref, linetype="dashed")+
        scale_x_continuous(limits=c(0, 300))+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)
}

plot_hcr_phase_diagram <- function(model_runs, extra_columns, dem_params, nyears){

    d <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params) %>%
        # SSB is females only
        filter(sex == "F", L1 == "naa") %>%
        # summarise SSB across year and sim 
        group_by(time, hcr, sim, L1) %>%
        summarise(spbio=sum(spbio)) %>%
        left_join(
            reshape2::melt(mse_small$hcr_f),
            by = c("time", "sim"),
        ) %>%
        group_by(time, hcr) %>%
        median_qi(spbio, value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        drop_na()

    ref_pts <- get_reference_points(model_runs, extra_columns, dem_params, year=nyears)

    ssbs <- seq(0, 300, 1)
    ssbs_df <- data.frame(spbio=ssbs, Fs=sapply(ssbs, \(x) npfmc_tier3_F(x, ref_pts$Bref, ref_pts$Fref)))

    plot <- ggplot(d)+
        geom_point(aes(x=spbio, y=value, color=time), size=1.5)+
        geom_segment(
            aes(x=spbio, y=value, xend=after_stat(lead(x)), yend=after_stat(lead(y)), color=time), 
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_line(data=ssbs_df, aes(x=spbio, y=Fs), size=0.75)+
        geom_hline(yintercept=ref_pts$Fref, linetype="dashed")+
        geom_vline(xintercept=ref_pts$Bref, linetype="dashed")+
        scale_x_continuous(limits=c(0, 300))+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)

}
