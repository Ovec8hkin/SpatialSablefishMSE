plot_ssb <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, base_hcr="F40"){
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

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- d %>% left_join(traj, by=traj_column) %>% filter(L1=="naa", hcr==hcr1) %>% group_by(om) %>% filter(time <= common)

    base_hcr_d <- d %>% filter(L1 == "naa", hcr == base_hcr)

    plot <- ggplot(d %>% filter(L1 == "naa")) + 
        geom_lineribbon(data = base_hcr_d, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        # geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        scale_y_continuous(limits=c(0, 500))+
        labs(x="Year", y="SSB")+
        coord_cartesian(expand=0)+
        guides(color=guide_legend(title="HCR", nrow=2), fill="none")+
        theme_bw()

    if(show_est){
        plot <- plot + geom_pointrange(data = d %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(.data[[v3]]))+guides(fill="none")
    }

    return(plot+custom_theme)
}

plot_relative_ssb <- function(data, v1="hcr", v2=NA, common_trajectory=64, base_hcr="No Fishing"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio")]
    
    base_ssb_data <- data %>% filter(hcr == base_hcr)
    rel_ssb <- data %>% left_join(base_ssb_data, by=c("time", "sim", "L1", "om"), suffix=c("", ".nofish")) %>%
        filter(time > common_trajectory) %>%
        mutate(
            rel_ssb = spbio/spbio.nofish
        ) %>%
        filter(L1 == "naa") %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(rel_ssb, .width=interval_widths)

    ylim <- c(0, 1.1)

    plot <- ggplot(rel_ssb) +
        geom_line(aes(x=time, y=rel_ssb, color=.data[[v1]], group=.data[[v1]]), size=0.85)+
        scale_y_continuous(limits=ylim)+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(expand=0)+
        labs(x="Year", y="Relative SSB")+
        guides(color=guide_legend(title="HCR", nrow=2))+
        facet_wrap(~.data[[v2]])

    return(plot+custom_theme)

}

plot_fishing_mortalities <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64){
    # Plot fishing mortality rates from OM and EM
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "F", "total_F")]

    f <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(F, total_F, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns)) %>%
        filter(name == "total_F")

    hcr1 <- as.character((f %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- f %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- f %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common)


    plot <- ggplot(f %>% filter(time > common_trajectory-1)) + 
        geom_lineribbon(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]))+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        guides(fill="none")+
        theme_bw()

    if(show_est){
        plot <- plot + geom_pointrange(data = f %>% filter(L1 == "faa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(.data[[v3]]))+guides(fill="none")
    }

    return(plot)
}

plot_recruitment <- function(data, v1="hcr", v2=NA, show_est=FALSE, common_trajectory=64){
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
        # geom_pointrange(data = r %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)+
        geom_hline(yintercept = mean_rec, linetype="dashed") + 
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        scale_y_continuous(limits=c(0, 120))+
        coord_cartesian(expand=0)+
        theme_bw()

    if(show_est){
        plot <- plot + geom_pointrange(data = f %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    if(!is.na(v2)){
        plot <- plot + facet_wrap(~.data[[v2]])
    }

    return(plot)
}

plot_landed_catch <- function(data, v1="hcr", v2=NA, v3=NA, by_fleet=FALSE, common_trajectory=64, base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(catch, total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((c %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- c %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    if(by_fleet){
        c <- c %>% filter(name == "catch")
    }else{
        c <- c %>% filter(name == "total_catch")
    }

    common <- c %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common)

    base_hcr_c <- c %>% filter(hcr == base_hcr)

    plot <- ggplot(c %>% left_join(traj, by=traj_column) %>% filter(time > common-1))+
        geom_lineribbon(data = base_hcr_c, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+ 
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        # scale_y_continuous(limits=c(0, 60))+
        labs(x="Year", y="Catch (mt)", color="HCR")+
        coord_cartesian(expand=0, ylim=c(0, 60))+
        guides(color=guide_legend(title="HCR", nrow=2), fill="none")+
        theme_bw()

    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(.data[[v3]]))+guides(fill="none")
    }
    
    # if(by_fleet){
    #     plot <- plot + facet_wrap(~fleet)
    # }

    return(plot+custom_theme)

}


plot_ssb_catch <- function(ssb_data, catch_data, v1="hcr", v2=NA, v3=NA, common_trajectory=64, base_hcr="F40"){

    group_columns <- colnames(ssb_data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio")]
    # Plot spawning biomass from OM and EM
    ssb_d <- ssb_data %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    hcr1 <- as.character((ssb_d %>% pull(hcr) %>% unique)[1])

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- ssb_d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    ssb_common <- ssb_d %>% left_join(traj, by=traj_column) %>% filter(L1=="naa", hcr==hcr1) %>% group_by(om) %>% filter(time <= common)


    group_columns <- colnames(catch_data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    catch_d <- catch_data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((catch_d %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- catch_d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    catch_common <- catch_d %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common)

    d <- bind_rows(ssb_d, catch_d) %>% filter(L1 != "naa_est") %>% 
            mutate(L1 = factor(L1, labels=c("Landed Catch", "SSB")))
    common <- bind_rows(ssb_common, catch_common) %>% filter(L1 != "naa_est") %>% 
            mutate(L1 = factor(L1, labels=c("Landed Catch", "SSB")))

    base_hcr <- d %>% filter(hcr == base_hcr)

    plot <- ggplot(d) + 
        geom_line(data = base_hcr, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        # geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        # scale_y_continuous(limits=c(0, 320))+
        labs(x="Year", y="1000s Metric Tons")+
        coord_cartesian(expand=0)+
        guides(color=guide_legend("HCR", nrow=2), fill="none")+
        theme_bw()+
        facet_grid(rows=vars(L1), cols=vars(.data[[v2]]), scales="free_y")+
        ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 500))
            )
        )
    return(plot+custom_theme)
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

plot_abc_tac <- function(data, v1="hcr", v2=NA, common_trajectory=64, base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "value")]

    q <- data %>%
        mutate(
            L1 = factor(L1, levels=c("abc", "tac", "exp_land", "attainment"), labels=c("ABC", "TAC", "Expected Landings", "Attainment"))
        ) %>%
        # filter(L1 != "Expected Landings") %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((q %>% pull(hcr) %>% unique)[1])
    traj_column <- v2
    traj <- q %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- q %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(eval(rlang::parse_expr(v2))) %>% filter(time <= common)

    base_hcr_q <- q %>% filter(hcr == base_hcr)

    plot <- ggplot(q %>% filter(time > common_trajectory-1))+
        geom_lineribbon(data = base_hcr_q, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=interaction(.data[[v1]], L1), color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(expand=0)+
        theme_bw()

    if(!is.na(v2)){
        plot <- plot + facet_grid(rows=vars(L1), cols=vars(.data[[v2]]), scales="free_y")
    }else{
        plot <- plot + facet_wrap(~L1, scales="free_y")
    }

    plot <- plot + ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 50)),
                scale_y_continuous(limits=c(0.5, 1.5))
            )
        )

    return(plot)
}

plot_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=64){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "total_F")]

    d <- data %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, total_F, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
            filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, total_F, hcr, om) %>% 
        rename(x=spbio, y=total_F) %>%
        group_by(hcr, om) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=total_F, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(data=ref_pts, aes(yintercept=Fref), linetype="dashed")+
        geom_vline(data=ref_pts, aes(xintercept=Bref), linetype="dashed")+
        scale_x_continuous(limits=c(0, 200))+
        scale_y_continuous(limits=c(0, 0.125))+
        coord_cartesian(expand=0)+
        facet_grid(cols=vars(hcr), rows=vars(om))

    return(plot)
}

plot_catch_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=64){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "total_catch")]

    d <- data %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, total_catch, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
            filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, total_catch, hcr, om) %>% 
        rename(x=spbio, y=total_catch) %>%
        group_by(hcr, om) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=total_catch, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(data=ref_pts, aes(yintercept=Fref), linetype="dashed")+
        geom_vline(data=ref_pts, aes(xintercept=Bref), linetype="dashed")+
        scale_x_continuous(limits=c(0, 200))+
        scale_y_continuous(limits=c(0, 35))+
        coord_cartesian(expand=0)+
        facet_grid(cols=vars(hcr), rows=vars(om))

    return(plot)
}

plot_hcr_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=64){

    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "value")]

    d <- data %>%
            rename(out_F=value) %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, out_F, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
            filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, out_F, hcr, om) %>% 
        rename(x=spbio, y=out_F) %>%
        group_by(hcr, om) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=out_F, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(data=ref_pts, aes(yintercept=Fref), linetype="dashed")+
        geom_vline(data=ref_pts, aes(xintercept=Bref), linetype="dashed")+
        scale_x_continuous(limits=c(0, 200))+
        scale_y_continuous(limits=c(0, 0.125))+
        coord_cartesian(expand=0)+
        facet_grid(cols=vars(hcr), rows=vars(om))

    return(plot)

}


plot_mse_summary <- function(model_runs, extra_columns, dem_params, hcr_filter, om_filter, common_trajectory=64){
    all_data <- bind_rows(
        get_ssb_biomass(model_runs, extra_columns, dem_params, hcr_filter, om_filter) %>% select(time, sim, L1, om, hcr, value=spbio),
        get_management_quantities(model_runs, extra_columns, hcr_filter, om_filter, spinup_years = common_trajectory),
        get_landed_catch(model_runs, extra_columns, hcr_filter, om_filter) %>% select(time, sim, L1, om, hcr, value=total_catch),
        get_fishing_mortalities(model_runs, extra_columns, hcr_filter, om_filter) %>% select(time, sim, L1, om, hcr, value=total_F)
    )

    ad <- all_data %>% 
        filter(L1 %in% c("tac", "faa", "naa", "land_caa")) %>%
        group_by(time, L1, om, hcr) %>%
        median_qi(value, .width = interval_widths) %>%
        reformat_ggdist_long(n=4) %>%
        mutate(
            L1 = factor(
                L1,
                levels = c("tac", "land_caa", "faa", "naa"),
                labels = c("TAC", "Catch", "Fishing Mortality", "Spawning Biomass")
            )
        )
    
    hcr1 <- as.character((ad %>% pull(hcr) %>% unique)[1])
    traj <- ad %>% distinct(om) %>% mutate(common=common_trajectory)

    common <- ad %>% left_join(traj, by="om") %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common)

    plot <- ggplot(ad %>% filter(time > common_trajectory-1), aes(x=time, y=median, color=hcr))+
        geom_line(size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85, color="black")+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        facet_grid(cols=vars(om), rows=vars(L1), scales="free_y")+
        facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 50), breaks=seq(0, 50, 10)),
                scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 15)),
                scale_y_continuous(limits=c(0, 0.2), breaks=seq(0, 0.20, 0.05)),
                scale_y_continuous(limits=c(0, 320, breaks=seq(0, 320, 50)))
            )
        )+
        scale_x_continuous(limits=c(0, ad %>% pull(time) %>% max))+
        scale_color_manual(values=hcr_colors)+
        labs(y="", x="Year", color="HCR")+
        coord_cartesian(expand=0)+
        theme_bw()

    return(plot)
}

plot_performance_metric_summary <- function(perf_data, v1="hcr", v2="om", is_relative=FALSE, summary_hcr="F40"){

    metric_minmax = perf_data %>% group_by(name) %>% summarise(min=min(lower), max=max(upper))
    axis_scalar <- c(0.9, 1.1)

    summary <- perf_data %>% filter(.width == 0.50, hcr == summary_hcr)
        # group_by(across(all_of(c(summary_var, "name")))) %>% 
        # summarise(value = mean(median))

    plot <- ggplot(perf_data)+
                geom_vline(data=summary, aes(xintercept = median), color="black")+
                scale_shape_discrete()+
                scale_color_manual(values=hcr_colors)+
                # facet_wrap(vars(name), scales="free_x")+
                labs(y="", x="", shape="OM", color="HCR")+
                coord_cartesian(expand=0)+
                guides(shape="none", color=guide_legend(nrow=1))+
                theme_bw()+
                theme(
                    plot.margin = margin(0.25, 1, 0.25, 0.25, "cm"),
                    panel.spacing.x = unit(5, "cm"),
                    plot.title = element_text(size=18),
                    legend.spacing.x = unit(1.5, "cm")
                )

    if(is.character(v2)){
        plot <- plot + 
                    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=.data[[v1]], color=.data[[v1]], shape=.data[[v2]]), point_size=3, position="dodge")+
                    facet_grid(rows=vars(.data[[v2]]), cols=vars(name), scales="free_x")
    }else{
        plot <- plot + 
                    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=.data[[v1]], color=.data[[v1]]), point_size=3, position="dodge")+
                    facet_wrap(~name, scales="free_x")
    }

    if(!is_relative){
        plot <- plot + 
                ggh4x::facetted_pos_scales(
                    x = list(
                        scale_x_continuous(limits=c(0, 55)),
                        scale_x_continuous(limits=c(0, 550), breaks=c(0, 150, 300, 450)),
                        scale_x_continuous(limits=c(0, 0.06), breaks=c(0, 0.02, 0.04, 0.06)),
                        scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.50, 1.0)),
                        # scale_x_continuous(limits=c(0, 5)),
                        # scale_x_continuous(limits=c(2, 15)),
                        scale_x_continuous(limits=c(5, 12)),
                        scale_x_continuous(limits=c(0, 5)),
                        scale_x_continuous(limits=c(0, 15))
                    
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

    return(plot+custom_theme)
}

plot_ssb_paginate <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, base_hcr="F40"){
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

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- d %>% left_join(traj, by=traj_column) %>% filter(L1=="naa", hcr==hcr1) %>% group_by(om) %>% filter(time <= common, time >= 40) %>% select(-hcr)

    base_hcr_d <- d %>% filter(L1 == "naa", hcr == base_hcr) %>% select(-hcr)

    ps <- lapply(unlist(unique(c(om_names))), function(o){

        d2 <- d %>% filter(om == o)
        base_hcr_d2 <- base_hcr_d %>% filter(om == o)

        ggplot(d2 %>% filter(L1 == "naa")) + 
            geom_lineribbon(data = base_hcr_d2, aes(x=time, y=median, ymin=lower, ymax=upper, group=1), color="black", size=0.85)+
            geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
            geom_line(data = common, aes(x=time, y=median), size=0.85)+
            geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
            # geom_hline(yintercept=121.4611, linetype="dashed")+
            scale_fill_brewer(palette="Blues")+
            scale_y_continuous(limits=c(0, max(d2 %>% pull(median))*1.2))+
            scale_x_continuous(limits=c(40, max(base_hcr_d2 %>% pull(time))))+
            facet_wrap(~ hcr, ncol=4, nrow=7)+
            labs(x="Year", y="SSB", title=o)+
            coord_cartesian(expand=0)+
            guides(fill="none", color="none")+
            theme_bw()+
            custom_theme+
            theme(
                plot.title = element_text(size=20)
            )
    })

    ggsave(
        filename = "~/Desktop/ssb_paginated.pdf", 
        plot = marrangeGrob(ps, nrow=1, ncol=1), 
        width = 8.5, height = 11
    )

    return(ps)
}

plot_catch_paginate <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c2 <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(catch, total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((c2 %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- c2 %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

        c2 <- c2 %>% filter(name == "total_catch")

    common <- c2 %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common, time >= 40) %>% select(-hcr)

    base_hcr_c <- c2 %>% filter(hcr == base_hcr) %>% select(-hcr)

    ps <- lapply(unlist(unique(c(om_names))), function(o){

        c3 <- c2 %>% filter(om == o)
        base_hcr_c2 <- base_hcr_c %>% filter(om == o)

        ggplot(c3 %>% filter(L1 == "land_caa")) + 
            geom_lineribbon(data = base_hcr_c2, aes(x=time, y=median, ymin=lower, ymax=upper, group=1), color="black", size=0.85)+
            geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
            geom_line(data = common, aes(x=time, y=median), size=0.85)+
            geom_vline(data=common, aes(xintercept=common_trajectory), linetype="dashed")+
            # geom_hline(yintercept=121.4611, linetype="dashed")+
            scale_fill_brewer(palette="Blues")+
            # scale_y_continuous(limits=c(0, max(c3 %>% pull(median))*1.2))+
            scale_x_continuous(limits=c(40, max(base_hcr_c2 %>% pull(time))))+
            facet_wrap(~ hcr, ncol=4, nrow=7)+
            labs(x="Year", y="Catch (mt)", title=o)+
            coord_cartesian(expand=0, ylim=c(0, 60))+
            guides(fill="none", color="none")+
            theme_bw()+
            custom_theme+
            theme(
                plot.title = element_text(size=20)
            )
    })

    ggsave(
        filename = "~/Desktop/catch_paginated.pdf", 
        plot = marrangeGrob(ps, nrow=1, ncol=1), 
        width = 8.5, height = 11
    )

    return(ps)
}

set_hcr_colors <- function(hcrs){
    hcr_colors <- hue_pal()(length(hcrs))
    hcr_colors[which(hcrs == "No Fishing")] <- "#000000"
    names(hcr_colors) <- hcrs
    return(hcr_colors)
}


custom_theme <- theme_bw()+theme(
    panel.spacing.y = unit(0.5, "cm"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=14),
    axis.text = element_text(size=14),
    strip.text = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position = "bottom"
)
