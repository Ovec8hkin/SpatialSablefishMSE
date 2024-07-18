plot_ssb <- function(data, v1="hcr", v2=NA, show_est=FALSE, common_trajectory=64){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio")]
    # Plot spawning biomass from OM and EM
    d <- data %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    hcr1 <- as.character((d %>% pull(hcr) %>% unique)[1])
    om1 <- as.character((d %>% pull(om) %>% unique)[1])

    plot <- ggplot(d %>% filter(L1 == "naa", time > common_trajectory-1)) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = d %>% filter(L1=="naa", hcr==hcr1, time <= common_trajectory), aes(x=time, y=median), size=0.85)+
        geom_vline(xintercept=common_trajectory, linetype="dashed")+
        geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 320))+
        coord_cartesian(expand=0)+
        guides(fill="none")+
        theme_bw()

    if(show_est){
        plot <- plot + geom_pointrange(data = d %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }

    return(plot)
}

plot_fishing_mortalities <- function(data, v1="hcr", v2=NA, show_est=FALSE, common_trajectory=64){
    # Plot fishing mortality rates from OM and EM
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "F", "total_F")]

    f <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(F, total_F, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns)) %>%
        filter(name == "total_F")

    hcr1 <- as.character((d %>% pull(hcr) %>% unique)[1])
    om1 <- as.character((d %>% pull(om) %>% unique)[1])

    plot <- ggplot(f %>% filter(L1 == "faa", time > common_trajectory-1)) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]))+
        geom_line(data = f %>% filter(L1=="faa", hcr==hcr1, time <= common_trajectory), aes(x=time, y=median), size=0.85)+
        geom_vline(xintercept=64, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        guides(fill="none")+
        theme_bw()

    if(show_est){
        plot <- plot + geom_pointrange(data = f %>% filter(L1 == "faa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
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

plot_landed_catch <- function(data, v1="hcr", v2=NA, by_fleet=FALSE, common_trajectory=64){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(catch, total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((c %>% pull(hcr) %>% unique)[1])
    om1 <- as.character((c %>% pull(om) %>% unique)[1])

    if(by_fleet){
        c <- c %>% filter(name == "catch")
    }else{
        c <- c %>% filter(name == "total_catch")
    }

    plot <- ggplot(c %>% filter(time > common_trajectory-1))+
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]))+
        geom_line(data = c %>% filter(hcr==hcr1, time <= common_trajectory), aes(x=time, y=median), size=0.85)+
        geom_vline(xintercept=64, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_y_continuous(limits=c(0, 60))+
        coord_cartesian(expand=0)+
        guides(fill="none")+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }
    
    # if(by_fleet){
    #     plot <- plot + facet_wrap(~fleet)
    # }

    return(plot)

}

plot_atage_trajectory_ternary <- function(data, segments, col_names){
    axis_names = names(data)[6:8]
    return(
        ggplot(data, aes(x=.data[[col_names[1]]], y=.data[[col_names[2]]], z=.data[[col_names[3]]], color=hcr))+
            coord_tern(Tlim=c(0, 1), Llim=c(0, 1), Rlim=c(0, 1))+
            geom_point()+
            geom_segment(
                data = segments, 
                aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend, group=hcr),
                arrow=arrow(length = unit(3, "mm"))
            )+
            scale_T_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_L_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_R_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            facet_grid(rows=vars(om), cols=vars(hcr))+
            labs(x=axis_names[1], y=axis_names[2], z=axis_names[3])+
            theme_bw()+
            theme(
                legend.position="bottom",
                tern.axis.arrow.show = TRUE,
                panel.spacing.x = unit(1, "cm"),
                panel.grid.minor = element_line(color="white")
            )
    
    )

}

plot_atage_density_ternary <- function(data, col_names){
    axis_names <- names(data)[6:8]
    return(
        ggplot(data, aes(x=.data[[col_names[1]]], y=.data[[col_names[2]]], z=.data[[col_names[3]]]))+
            coord_tern()+
            stat_density_tern(
                geom='polygon',
                aes(fill=..level..),
                bins=100,
            )+
            geom_mean_ellipse(color="white")+
            scale_fill_viridis(limits=c(0, 65))+
            scale_T_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_L_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_R_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            facet_grid(rows=vars(om), cols=vars(hcr))+
            labs(fill="Simulation Years", x=axis_names[1], y=axis_names[2], z=axis_names[3])+
            theme_bw()+
            theme(
                legend.position="bottom",
                tern.axis.arrow.show = TRUE,
                panel.spacing.x = unit(1, "cm"),
                panel.grid.minor = element_line(color="white")
            )
    )
}

# plot_atage <- function(data, v1="hcr", v2=NA){
#     group_columns <- colnames(data)
#     group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]
    
#     d <- data %>%
#         group_by(time, class, hcr, om, L1) %>%
#         median_qi(value, .width=c(0.50))

#     plot <- ggplot()+
#         geom_bar(aes(x=time, y=value, fill=class), position="fill", stat="identity")+
#         geom_vline(xintercept = 2022-1960+1, color="white", size=1.25, linetype="dashed")+
#         scale_fill_viridis(direction=-1, discrete=TRUE, option="magma")+
#         scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
#         coord_cartesian(expand=0)+
#         labs(x="Year", fill="Age Group")+
#         guides(fill=guide_legend(reverse=TRUE))+
#         # facet_wrap(~L1, ncol=1, scales="free_x")+
#         facet_grid(rows=vars(hcr), cols=vars(L1))+
#         theme_bw()+
#         theme(
#             axis.text = element_text(size=12),
#             axis.title.y=element_blank(), 
#             strip.background = element_blank(),
#             strip.text.x = element_text(size=16, hjust=0),
#             panel.spacing.y = unit(0.4, "in"),
#             legend.position = "bottom"
#         )

#     if(!is.na(v2)){
#         plot <- plot + facet_wrap(~.data[[v2]])
#     }else{
    
#     }
    
#     # if(by_fleet){
#     #     plot <- plot + facet_wrap(~fleet)
#     # }

#     return(plot)

# }

plot_abc_tac <- function(data, v1="hcr", v2=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "value")]

    q <- data %>%
        # filter(L1 != "Expected Landings") %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    plot <- ggplot(q %>% filter(time > 63))+
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, color=.data[[v1]], group=interaction(.data[[v1]], L1)))+
        geom_line(data=q %>% filter(hcr=="F40", time <= 64), aes(x=time, y=median), color="black", size=0.85)+
        geom_vline(xintercept=64, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        coord_cartesian(expand=0)+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_grid(rows=vars(L1), cols=vars(.data[[v2]]), scales="free_y")
    }else{
        plot <- plot + facet_wrap(~L1, scales="free_y")
    }

    plot <- plot + ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 100)),
                scale_y_continuous(limits=c(0, 100)),
                scale_y_continuous(limits=c(0, 100)),
                scale_y_continuous(limits=c(0.5, 1))
            )
        )

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


plot_mse_summary <- function(model_runs, extra_columns, dem_params, common_trajectory=64){
    all_data <- bind_rows(
        get_ssb_biomass(model_runs, extra_columns, dem_params) %>% select(time, sim, L1, om, hcr, value=spbio),
        get_management_quantities(model_runs, extra_columns),
        get_landed_catch(model_runs, extra_columns) %>% select(time, sim, L1, om, hcr, value=total_catch),
        get_fishing_mortalities(model_runs, extra_columns) %>% select(time, sim, L1, om, hcr, value=total_F)
    )

    ad <- all_data %>% 
        filter(L1 %in% c("abc", "faa", "naa", "land_caa")) %>%
        group_by(time, L1, om, hcr) %>%
        median_qi(value, .width = interval_widths) %>%
        reformat_ggdist_long(n=4) %>%
        mutate(
            L1 = factor(
                L1,
                levels = c("abc", "land_caa", "faa", "naa"),
                labels = c("ABC", "Catch", "Fishing Mortality", "Spawning Biomass")
            )
        )

    plot <- ggplot(ad %>% filter(time > common_trajectory-1), aes(x=time, y=median, color=hcr))+
        geom_line(size=0.85)+
        geom_vline(xintercept = common_trajectory, linetype="dashed")+
        geom_line(data=ad %>% filter(hcr=="F40", time <= common_trajectory), color="black", size=0.85)+
        facet_grid(cols=vars(om), rows=vars(L1), scales="free_y")+
        facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, 25)),
                scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 15)),
                scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.20, 0.05)),
                scale_y_continuous(limits=c(0, 650, breaks=seq(0, 650, 100)))
            )
        )+
        scale_x_continuous(limits=c(0, ad %>% pull(time) %>% max))+
        labs(y="", x="Year", color="HCR")+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)
}

plot_performance_metric_summary <- function(perf_data, is_relative=FALSE){

    om_summary <- perf_data %>% filter(.width == 0.50) %>% 
        group_by(om, name) %>% 
        summarise(value = mean(median))

    plot <- ggplot(perf_data)+
                geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=hcr, color=hcr, shape=om), point_size=3, position="dodge")+
                geom_vline(data=om_summary, aes(xintercept = value), color="black")+
                scale_shape_discrete()+
                facet_grid(rows=vars(om), cols=vars(name), scales="free_x")+
                # facet_wrap(vars(name), scales="free_x")+
                labs(y="", x="", shape="OM", color="HCR")+
                coord_cartesian(expand=0)+
                guides(shape=guide_legend(nrow=3), color=guide_legend(nrow=4))+
                theme_bw()+
                theme(
                    plot.margin = margin(0.25, 1, 0.25, 0.25, "cm"),
                    panel.spacing.x = unit(0.5, "cm"),
                    plot.title = element_text(size=18),
                    legend.spacing.x = unit(1.5, "cm")
                )

    if(!is_relative){
        plot <- plot + 
            ggh4x::facetted_pos_scales(
                x = list(
                    scale_x_continuous(limits=c(0, 65), breaks=seq(0, 60, 20), labels = seq(0, 60, 20)),
                    scale_x_continuous(limits=c(0, 600), breaks=seq(0, 600, 200), labels=seq(0, 600, 200)),
                    scale_x_continuous(limits=c(0, 0.07), breaks=seq(0, 0.06, 0.02), labels=seq(0, 0.06, 0.02)),
                    scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25), labels=seq(0, 100, 25)),
                    scale_x_continuous(limits=c(0, 0.25), breaks=seq(0, 0.25, 0.05), labels=seq(0, 25, 5)),
                    scale_x_continuous(limits=c(0, 45), breaks=seq(0, 40, 10), labels=seq(0, 40, 10))
                )
            )
    }else{
        plot <- plot + ggh4x::facetted_pos_scales(
                x = list(
                    scale_x_continuous(limits=c(0, 1.5)),
                    scale_x_continuous(limits=c(0.75, 2.0)),
                    scale_x_continuous(limits=c(0, 3.5)),
                    scale_x_continuous(limits=c(0.75, 1.5)),
                    scale_x_continuous(limits=c(0.75, 2.5)),
                    scale_x_continuous(limits=c(0.5, 1.25))
                )
            )
    }

    return(plot)
}


custom_theme <- theme(
    panel.spacing.y = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position = "bottom"
)
