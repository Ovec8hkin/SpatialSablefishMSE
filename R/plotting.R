plot_ssb <- function(data, v1="hcr", v2=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio")]
    # Plot spawning biomass from OM and EM
    d <- data %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))


    plot <- ggplot(d %>% filter(L1 == "naa")) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.4)+
        geom_pointrange(data = d %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_vline(xintercept=64, linetype="dashed")+
        geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 320))+
        coord_cartesian(expand=0)+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])
    }

    return(plot)
}

plot_fishing_mortalities <- function(data, v1="hcr", v2=NA){
    # Plot fishing mortality rates from OM and EM
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "F", "total_F")]

    f <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(F, total_F, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns)) %>%
        filter(name == "total_F")

    plot <- ggplot(f %>% filter(L1 == "faa")) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.4)+
        geom_pointrange(data = f %>% filter(L1 == "faa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_vline(xintercept=64, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])
    }

    return(plot)
}

plot_recruitment <- function(data, v1="hcr", v2=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "rec")]

    r <- data %>%
        # summarise SSB across year and sim 
        group_by(across(all_of(group_columns))) %>%
        median_qi(rec, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))

    mean_rec <- r %>% pull(median) %>% mean

    plot <- ggplot(r %>% filter(L1 == "naa")) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.4)+
        geom_pointrange(data = r %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_hline(yintercept = mean_rec, linetype="dashed") + 
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 120))+
        coord_cartesian(expand=0)+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])
    }

    return(plot)
}

plot_landed_catch <- function(data, v1="hcr", v2=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(catch, total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    if(by_fleet){
        c <- c %>% filter(name == "catch")
    }else{
        c <- c %>% filter(name == "total_catch")
    }

    plot <- ggplot(c, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]))+
        geom_lineribbon()+
        geom_vline(xintercept=64, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 50))+
        coord_cartesian(expand=0)+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])
    }
    
    # if(by_fleet){
    #     plot <- plot + facet_wrap(~fleet)
    # }

    return(plot)

}

plot_abc_tac <- function(data, v1="hcr", v2=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "value")]

    q <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    plot <- ggplot(q)+
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, color=.data[[v1]], group=interaction(.data[[v1]], L1)))+
        scale_y_continuous(limits=c(0, 70))+
        scale_fill_brewer(palette="Blues")+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_grid(rows=vars(L1), cols=vars(.data[[v2]]))
    }else{
        plot <- plot + facet_wrap(~L1)
    }

    return(plot)
}

plot_phase_diagram <- function(model_runs, extra_columns, dem_params, nyears){
    d <- get_ssb_biomass(model_runs, extra_columns) %>%
        # SSB is females only
        filter(sex == "F", L1 == "naa") %>%
        # summarise SSB across year and sim 
        group_by(time, hcr, sim, L1) %>%
        summarise(spbio=sum(spbio)) %>%
        left_join(
            bind_mse_outputs(model_runs, "out_f", extra_columns),
            by = c("time", "sim", "hcr"),
        ) %>%
        group_by(time, hcr) %>%
        median_qi(spbio, value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        filter(.width == 0.50)

    ref_pts <- get_reference_points(model_runs[1], list(hcr="test"), dem_params, year=nyears)

    segments <- d %>% as_tibble() %>% 
        select(spbio, value, hcr) %>% 
        rename(x=spbio, y=value) %>%
        group_by(hcr) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=value, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
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
            bind_mse_outputs(model_runs, "hcr_f", extra_columns),
            by = c("time", "sim", "hcr"),
        ) %>%
        group_by(time, hcr) %>%
        median_qi(spbio, value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        drop_na() %>%
        filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, value, hcr) %>% 
        rename(x=spbio, y=value) %>%
        group_by(hcr) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    ref_pts <- get_reference_points(model_runs[1], list(hcr="test"), dem_params, year=nyears)

    ssbs <- seq(0, 300, 1)
    ssbs_df <- data.frame(spbio=ssbs, Fs=sapply(ssbs, \(x) npfmc_tier3_F(x, ref_pts$Bref, ref_pts$Fref)))

    plot <- ggplot(d)+
        geom_point(aes(x=spbio, y=value, color=hcr, group=hcr), size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, color=hcr, group=hcr),
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
