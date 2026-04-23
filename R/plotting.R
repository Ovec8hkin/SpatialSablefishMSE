plot_ssb <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=54, time_horizon=NULL, interval_widths=c(0.50, 0.80), base_hcr="F40", relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "biomass")]
    d <- data %>% distinct(.keep_all=TRUE) %>% select(-biomass)

    if(!is.null(time_horizon)){
        d <- d %>% filter_times(time_horizon)
    }

    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        d <- d %>% 
            group_by(across(all_of(c(groups, "sim")))) %>% 
            summarise(
                spbio=sum(spbio)
            ) %>% 
            mutate(region="Alaska")
        v3 <- "region"
    }

    # Relativize to specific HCR
    hcrs <- d %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        d <- d %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="spbio", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    # Plot spawning biomass from OM and EM
    d <- d %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=interval_widths, .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))

    plot <- plot_timeseries(d %>% filter(L1 == "naa"), v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="SSB (1000s mt)")
    
    if(show_est){
        plot <- plot + geom_pointrange(
                            data = d %>% filter(L1 == "naa_est"), 
                            aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), 
                            alpha=0.35
                        )
    }
    
    if(!is.na(relative)){
        plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    }

    return(plot)
}

plot_average_age <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, time_horizon=NULL, interval_widths=c(0.50, 0.80), base_hcr="F40", relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "avg_age")]
    d <- data %>% distinct(.keep_all=TRUE)

    if(!is.null(time_horizon)){
        d <- d %>% filter_times(time_horizon)
    }

    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        d <- d %>% group_by(across(all_of(c(groups, "sim")))) %>% summarise(avg_age = mean(avg_age)) %>% mutate(region="Alaska")
        v3 <- "region"
    }

    # Relativize to specific HCR
    hcrs <- d %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        d <- d %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="avg_age", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    # Plot spawning biomass from OM and EM
    d <- d %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(avg_age, .width=interval_widths, .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    plot <- plot_timeseries(d, v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Average Age (Years)")
    if(!is.na(relative)){
        plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    }

    return(plot)
}

plot_numbers_at_pricegrade <- function(data, v1="hcr", v2=NA, v3=NA, price_grade="Grade 1/2 (1-2yo)", common_trajectory=54, time_horizon=NULL, interval_widths=c(0.50, 0.80), base_hcr="F40", relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "value")]
    d <- data %>% distinct(.keep_all=TRUE)

    if(!is.null(time_horizon)){
        d <- d %>% filter_times(time_horizon)
    }

    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        d <- d %>% filter(class == price_grade) %>% 
            group_by(across(all_of(c(groups, "sim")))) %>% 
            summarise(value = sum(value)) %>% 
            mutate(region="Alaska")
        v3 <- "region"
    }

    # Relativize to specific HCR
    hcrs <- d %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        d <- d %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="value", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    # Plot spawning biomass from OM and EM
    d <- d %>%
        filter(class == price_grade) %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(value, .width=interval_widths, .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    plot <- plot_timeseries(d, v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Average Age (Years)")
    if(!is.na(relative)){
        plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    }

    return(plot)
}

plot_depletion <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40", scales="fixed"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "biomass")]
    d <- data
    if(is.na(v3)){
        group_columns <- group_columns[group_columns != "region"]
        d <- d %>% group_by(across(all_of(c(group_columns, "sim")))) %>% summarise(spbio=sum(spbio), biomass=sum(biomass))
    }
    # Plot spawning biomass from OM and EM
    d <- d %>%
        select(-c("biomass")) %>%
        group_by(across(all_of(group_columns[group_columns != "time"]))) %>%
        mutate(dep=spbio/spbio[time==1]) %>%
        ungroup() %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(dep, .width=interval_widths, .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    return(
        plot_timeseries(d, v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Average Age (Years)")
    )
}

plot_relative_ssb <- function(data, v1="hcr", v2=NA, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="No Fishing"){
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
        guides(color=guide_legend(title="Management \n Strategy", nrow=2))+
        facet_wrap(~.data[[v2]])

    return(plot+custom_theme)

}

plot_fishing_mortalities <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=54, interval_widths=c(0.50, 0.80)){
    # Plot fishing mortality rates from OM and EM
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "F")]

    f <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(F, .width=interval_widths, .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))

    hcr1 <- as.character((f %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- f %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- f %>% left_join(traj, by=traj_column) %>% 
                filter(hcr==hcr1, time <= common, .width == 0.50, L1=="faa")

    base_hcr_f <- f %>% filter(L1 == "faa", hcr == hcr1)

    plot <- ggplot(f %>% filter(time > common_trajectory-1, L1 == "faa")) + 
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=interaction(.data[[v1]], fleet), linetype=fleet, color=.data[[v1]]))+
        geom_line(data = common, aes(x=time, y=median, linetype=fleet), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=c("Red", "Blue"))+
        scale_y_continuous(limits=c(0, 0.20))+
        coord_cartesian(expand=0)+
        guides(fill="none")

    if(show_est){
        plot <- plot + geom_pointrange(data = f %>% filter(L1 == "faa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(.data[[v3]]))+guides(fill="none")
    }

    return(plot+custom_theme)
}

plot_recruitment <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=54, interval_widths=c(0.50, 0.80), time_horizon=NULL, base_hcr="F40", relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "rec")]
    r <- data %>% distinct(.keep_all=TRUE)

    if(!is.null(time_horizon)){
        r <- r %>% filter_times(time_horizon)
    }

    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        r <- r %>% 
            group_by(across(all_of(c(groups, "sim")))) %>% 
            summarise(
                rec=sum(rec)
            ) %>% 
            mutate(region="Alaska")
        v3 <- "region"

        mean_rec <- r %>% filter(time > common_trajectory, L1=="naa") %>%
            group_by(om) %>%
            summarise(
                rec = mean(rec, na.rm=TRUE)
            )
    }else{
        mean_rec <- r %>% filter(time > common_trajectory, L1=="naa") %>%
            group_by(om, region) %>%
            summarise(
                rec = mean(rec, na.rm=TRUE)
            )
    }

    # Relativize to specific HCR
    hcrs <- r %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        r <- r %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="rec", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    r <- r %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(rec, .width=interval_widths, .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))

    max_year <- r %>% pull(time) %>% max

    plot <- plot_timeseries(r %>% filter(L1 == "naa"), v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Recruitment (millions)")
    
    if(show_est){
        plot <- plot + geom_pointrange(
                            data = r %>% filter(L1 == "naa_est"), 
                            aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), 
                            alpha=0.35
                        )
    }

    return(plot)
}

plot_random_recruitment_trajectories <- function(data, seed_list, nsims, time_horizon, by_region=FALSE){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("rec")]
    
    groups <- group_columns
    if(!by_region){
        groups <- group_columns[group_columns != "region"]
    }

    r <- recruit_data %>%
        filter_times(time_horizon) %>%
        filter(sim %in% sample(seed_list, nsims)) %>%
        group_by(across(all_of(groups))) %>%
        summarise(rec = sum(rec, na.rm=TRUE))

    ymax <- max(r$rec)

    mean_rec <- r %>% group_by(om) %>% summarise(r=median(rec))

    summ_rec <- r %>%
        group_by(time, om) %>%
        median_qi(rec, .width=interval_widths)

    plot <- ggplot(r %>% filter(L1 == "naa"))+
        geom_line(aes(x=time, y=rec, color=hcr, fill=hcr, group=sim), alpha=0.8)+
        geom_line(data=summ_rec, aes(x=time, y=rec, ymin=.lower, ymax=.upper), color="black", size=1, alpha=0.75)+
        # geom_line(aes(x=time, y=rec, color=om, group=sim), size=0.5, alpha=0.6)+
        geom_hline(data=mean_rec, aes(yintercept=r), linetype="dashed")+
        geom_text(data=mean_rec, aes(x=100, y=ymax*0.90, label=round(r, 3)), size=6)+
        labs(x="Year", y="Recruitment (millions)", color="Management \nStrategy")+
        custom_theme

    if(!by_region){
        plot <- plot + facet_grid(hcr~om)
    }else{
        plot <- plot + facet_grid(rows=vars(region), cols=vars(om), scales="free_y")
    }

    return(plot)
}

plot_landed_catch <- function(data, v1="hcr", v2=NA, v3=NA, by_fleet=FALSE, common_trajectory=54, time_horizon=NULL, interval_widths=c(0.50, 0.80), base_hcr="F40", relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]
    c <- data

    if(!is.null(time_horizon)){
        c <- c %>% filter_times(time_horizon)
    }

    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        c <- c %>% 
            group_by(across(all_of(c(groups, "sim")))) %>% 
            summarise(
                catch=sum(catch), 
                total_catch=sum(total_catch)
            ) %>% 
            mutate(region="Alaska")
        v3 <- "region"
    }

    val <- ifelse(by_fleet, "catch", "total_catch")
    c <- c %>% select(all_of(c(group_columns, "sim", val)))

    # Relativize to specific HCR
    hcrs <- c %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        c <- c %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column=val, 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    c <- c %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(eval(rlang::parse_expr(val)), .width=interval_widths, .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))


    plot <- plot_timeseries(c, v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Catch (1000s mt)")
    # if(!is.na(relative)){
    #     plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    # }

    return(plot)

}

plot_econ_value <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, time_horizon=NULL, interval_widths=c(0.50, 0.80), base_hcr="F40", relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "total_value")]
    d <- data %>% distinct(.keep_all=TRUE)

    if(!is.null(time_horizon)){
        d <- d %>% filter_times(time_horizon)
    }

    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        d <- d %>% group_by(across(all_of(c(groups, "sim")))) %>% summarise(total_value=sum(total_value)) %>% mutate(region="Alaska")
        v3 <- "region"
    }

    # Relativize to specific HCR
    hcrs <- d %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        d <- d %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="total_value", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    # Plot spawning biomass from OM and EM
    d <- d %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(total_value, .width=interval_widths, .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    plot <- plot_timeseries(d, v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Economic Value")
    if(!is.na(relative)){
        plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    }

    return(plot)
}

plot_ssb_catch <- function(ssb_data, catch_data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40"){

    group_columns <- colnames(ssb_data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "biomass")]
    # Plot spawning biomass from OM and EM
    ssb_d <- ssb_data %>%
        select(-c("biomass")) %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=interval_widths, .simple_names=FALSE) %>%
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
        median_qi(total_catch, .width=interval_widths, .simple_names=TRUE) %>%
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
        geom_line(data = base_hcr, aes(x=time, y=median, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        # geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        # scale_y_continuous(limits=c(0, 320))+
        labs(x="Year", y="1000s Metric Tons")+
        coord_cartesian(expand=0)+
        guides(color=guide_legend("Management \n Strategy", nrow=2), fill="none")+
        facet_grid(rows=vars(L1), cols=vars(.data[[v2]]), scales="free_y")+
        ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 500))
            )
        )
    return(plot+custom_theme)
}

plot_abc_tac <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40", by_fleet=FALSE){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "value")]

    q <- data
    if(is.na(v3)){
        groups <- group_columns[group_columns != "region"]
        q <- q %>% group_by(across(all_of(c(groups, "sim")))) %>% summarise(value=sum(value)) %>% mutate(region="Alaska")
        v3 <- "region"
    }

    if(!by_fleet){
        q <- q %>%
            group_by(across(all_of(c(group_columns[group_columns != "fleet"], "sim")))) %>% 
            summarise(value=sum(value))
        group_columns <- group_columns[! group_columns %in% c("fleet")]
    }

    q <- q %>%
        mutate(
            L1 = factor(L1, levels=c("abc", "tac", "exp_land"), labels=c("ABC", "TAC", "Expected Landings"))
        ) %>%
        # filter(L1 != "Expected Landings") %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(value, .width=interval_widths, .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((q %>% pull(hcr) %>% unique)[1])
    traj_column <- v2
    traj <- q %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- q %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(eval(rlang::parse_expr(v2))) %>% filter(time <= common)

    base_hcr_q <- q %>% filter(hcr == base_hcr)

    plot <- ggplot(q)+
        # geom_lineribbon(data = base_hcr_q, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, color=.data[[v1]], linetype=L1, group=.data[[v1]]), size=0.85)+
        # geom_line(data = common, aes(x=time, y=median, linetype=L1), size=0.85)+
        # geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(expand=0)+
        guides(fill="none")


    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + ggh4x::facet_nested(
            rows=vars(.data[[v2]], fleet), 
            cols=vars(.data[[v3]]), 
            scales="free_y"
        )
    }

    return(plot+custom_theme)
}

#' Plot Estimation Bias
#' Create either a boxplot or a timeseries plot of estimation bias across
#' OMs and HCRs.
#' 
#' @param bias_data tibble produced from get_estimation_bias()
#' @param type either "boxplot" or "timeseries"
#' 
#' @export plot_estimation_bias
#' 
plot_estimation_bias <- function(bias_data, type="boxplot"){
    if(type == "boxplot"){
        plot <- ggplot(bias_data)+
            geom_boxplot(aes(x=om, y=bias))+
            geom_hline(yintercept=0, linetype="dashed")+
            facet_wrap(~hcr)+
            scale_x_discrete(labels=scales::label_wrap(15))+
            custom_theme
    }else if(type=="timeseries"){
        plot <- ggplot(data=bias_data, aes(x=time, y=bias, color=hcr, fill=hcr))+
            stat_lineribbon(.width=c(0.50, 0.80), alpha=0.25)+
            geom_hline(yintercept=0, linetype="dashed")+
            facet_wrap(~om)+
            custom_theme
    }
    return(plot)
}

plot_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=54, interval_widths=c(0.50, 0.80)){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "total_F")]

    d <- data %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, total_F, .width=interval_widths, .simple_names=FALSE) %>%
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

plot_catch_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=54, interval_widths=c(0.50, 0.80)){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "total_catch")]

    d <- data %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, total_catch, .width=interval_widths, .simple_names=FALSE) %>%
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

plot_hcr_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=54, interval_widths=c(0.50, 0.80)){

    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "value")]

    d <- data %>%
            rename(out_F=value) %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, out_F, .width=interval_widths, .simple_names=FALSE) %>%
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

plot_performance_metric_summary <- function(perf_data, v1="hcr", v2="om", is_relative=FALSE, summary_hcr="F40"){

    metric_minmax = perf_data %>% group_by(name) %>% summarise(min=min(lower), max=max(upper))
    axis_scalar <- c(0.9, 1.1)

    summary <- perf_data %>% filter(!is.infinite(median), hcr != "No Fishing") %>% summarise(median=mean(median))

    plot <- ggplot(perf_data)+
                geom_vline(data=summary, aes(xintercept = median), color="black")+
                scale_shape_discrete()+
                scale_color_manual(values=rank_colors)+
                # scale_color_manual(values=hcr_colors)+
                # facet_wrap(vars(name), scales="free_x")+
                labs(y="", x="", shape="OM", color="Performance Order")+
                coord_cartesian(expand=0)+
                guides(shape="none", color=guide_legend(nrow=1))
                theme(
                    plot.margin = margin(0.25, 1, 0.25, 0.25, "cm"),
                    panel.spacing.x = unit(5, "cm"),
                    plot.title = element_text(size=18),
                    legend.spacing.x = unit(1.5, "cm")
                )

    if(is.character(v2)){
        plot <- plot + 
                    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=.data[[v1]], color=rank, shape=.data[[v2]]), point_size=3, position="dodge")+
                    facet_grid(rows=vars(.data[[v2]]), cols=vars(name), scales="free_x")
    }else{
        plot <- plot + 
                    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=.data[[v1]], color=rank), point_size=3, position="dodge")+
                    facet_wrap(~name, scales="free_x")
    }

    if(!is_relative){
        plot <- plot + 
                ggh4x::facetted_pos_scales(
                    x = list(
                        scale_x_continuous(limits=c(0, 55)),
                        scale_x_continuous(limits=c(0, 0.06), breaks=c(0, 0.02, 0.04, 0.06)),
                        # scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.50, 1.0)),
                        scale_x_continuous(limits=c(0, 550), breaks=c(0, 150, 300, 450)),
                        scale_x_continuous(limits=c(0, 15)),
                        scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2))
                    )
                )
    }else{
        plot <- plot + ggh4x::facetted_pos_scales(
                x = list(
                    scale_x_continuous(limits=c(0, 2)),
                    scale_x_continuous(limits=c(0, 1.25)),
                    scale_x_continuous(limits=c(0, 3.5)),
                    scale_x_continuous(limits=c(0, 2)),
                    scale_x_continuous(limits=c(0.75, 2.5)),
                    scale_x_continuous(limits=c(0.5, 1.25))
                )
            )
    }

    return(plot+custom_theme)
}

plot_ssb_paginate <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio")]
    # Plot spawning biomass from OM and EM
    d <- data %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=interval_widths, .simple_names=FALSE) %>%
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
            guides(fill="none", color="none")
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

plot_catch_paginate <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c2 <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(catch, total_catch, .width=interval_widths, .simple_names=TRUE) %>%
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
            guides(fill="none", color="none")
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

#' Plot EM Diagnostics
#' 
#' Take model runs and plot set of model dignostics, including:
#' 1) Model fits to catch and indices
#' 2) Model fits to age composition data
#' 3) Model fits to derived quantities (SSB, harvest rate, recruitment)
#'
#' @param hcr name of HCR to plot diagnostics for
#' @param om name of OM to plot diagnostics for
#' @param simulation_seed simualtion seed to plot diagnostics for
#' @param simulation_year Year of the simulation to plot diagnostics for
#' @param report_convergence whether to calculate and report simulation convergence
#' 
#' @return ggplot figure showing model diagnostics plots
#' 
#' @export plot_diags
#' 
plot_diags <- function(hcr, om, simulation_seed, simulation_year, om_list, seed_list, report_convergence=TRUE){
    
    simulation_year <- simulation_year-1 # simulation years offset by 1

    x <- find_model_run(hcr, om, simulation_seed, om_list, seed_list)
    m <- x$m
    mse_obj <- x$model_run
    
    n_proj_years <- m$mse_options_list$mse_options$n_proj_years
    n_spinup_years <- m$mse_options_list$mse_options$n_spinup_years
    seeds <- m$seeds
    nsims <- length(seeds)
    simulation_number <- which(seeds == simulation_seed)
    
    em_model_obj <- mse_obj$model_outs[[(simulation_year+1)+(simulation_number-1)*(n_proj_years+1)]]

    report <- em_model_obj$rep

    if(report_convergence){
        convergence_message <- capture.output(
            tryCatch({
                sdrep = TMB::sdreport(em_model_obj)
                SPoRC::post_optim_sanity_checks(sdrep, report)
            }, error = function(e){
                # cat(e$message)
            }),
            type="message"
        )
        split <- str_split(convergence_message, "\\. ")[[1]]
        convergence_message <- split[-length(split)]
    }

    dem_params <- mse_obj$om$dem_params

    nyears <- n_spinup_years+simulation_year
    naa <- mse_obj$naa[1:nyears,,,,simulation_number]
    rec <- mse_obj$naa[1:nyears,1,,,simulation_number]
    faa <- mse_obj$faa[1:nyears,,,,,simulation_number]
    catch <- mse_obj$land_caa[1:nyears,,,,,simulation_number]
    survey_obs <- mse_obj$survey_obs

    # Plot SSB Fit
    em_ssb <- as.vector(report$SSB)
    om_ssb <- apply(naa[,,1,,drop=FALSE]*dem_params$mat[1:nyears,,1,,drop=FALSE]*dem_params$waa[1:nyears,,1,,drop=FALSE], 1, sum)

    ssb_df <- tibble(year=1:nyears, om=om_ssb, em=em_ssb)

    ssb_plot <- ggplot(ssb_df)+
        geom_line(aes(x=year, y=em))+
        geom_line(aes(x=year, y=om), color="red", alpha=0.30)+
        geom_point(aes(x=year, y=om), color="red")+
        labs(x="Simulation Year", y="SSB", title="SSB Comparison")+
        custom_theme

    # Plot harvest rate fit
    em_catch <- apply(report$PredCatch, 2, sum)
    em_prop_fs <- apply(aperm(report$FAA, c(2, 3, 4, 1, 5))[1:nyears,,,1,,drop=FALSE], c(1, 5), max)/apply(apply(aperm(report$FAA, c(2, 3, 4, 1, 5))[1:nyears,,,1,,drop=FALSE], c(1, 5), max), 1, sum)
    em_jointselret <- calculate_joint_selret(
        sel = aperm(report$fish_sel, c(2, 3, 4, 1, 5)),
        ret = dem_params$ret[1:nyears,,,1,,drop=FALSE],
        prop_fs = em_prop_fs
    )
    em_exploit_bio <- apply(
        aperm(report$NAA, c(2, 3, 4, 1))[1:nyears,,,,drop=FALSE]*em_jointselret$sel*dem_params$waa[1:nyears,,,1,drop=FALSE],
        1,
        sum
    )
    em_harvest_rate <- em_catch/em_exploit_bio


    om_catch <- apply(catch, 1, sum)
    om_prop_fs <- apply(faa[1:nyears,,,,,drop=FALSE], c(1, 4, 5), max)/apply(apply(faa[1:nyears,,,,,drop=FALSE], c(1, 4, 5), max), 1, sum)
    om_jointselret <- calculate_joint_selret(
        sel = dem_params$sel[1:nyears,,,,,drop=FALSE],
        ret = dem_params$ret[1:nyears,,,,,drop=FALSE],
        prop_fs = om_prop_fs
    )
    om_exploit_bio <- apply(
        naa[1:nyears,,,,drop=FALSE]*om_jointselret$sel*dem_params$waa[1:nyears,,,,drop=FALSE],
        1,
        sum
    )
    om_harvest_rate <- om_catch/om_exploit_bio

    harvest_rate_df <- data.frame("time"=1:nyears, om=om_harvest_rate, em=em_harvest_rate)
    hr_plot <- ggplot(harvest_rate_df)+
        geom_line(aes(x=time, y=em))+
        geom_line(aes(x=time, y=om), color="red", alpha=0.30)+
        geom_point(aes(x=time, y=om), color="red")+
        labs(x="Simulation Year", y="Harvest Rate", title="Harvest Rate Comparison")+
        custom_theme


    # Plot recruitment fit
    em_rec <- as.vector(report$Rec)
    om_rec <- apply(naa[1:nyears,1,,,drop=FALSE], 1, sum)

    rec_df <- tibble(year=1:nyears, om=om_rec, em=em_rec)
    rec_plot <- ggplot(rec_df)+
        geom_line(aes(x=year, y=em))+
        geom_line(aes(x=year, y=om), color="red", alpha=0.30)+
        geom_point(aes(x=year, y=om), color="red")+
        labs(x="Simulation Year", y="Recruits", title="Recruitment Comparison")+
        custom_theme

    # Plot survey fits
    em_survey <- aperm(report$PredSrvIdx, c(2, 1, 3))[1:nyears,,1:2]
    dimnames(em_survey) <- list("time"=1:nyears, "fleet"=c("LLS", "BTS"))

    om_ll_survey <- apply(survey_obs$rpns[1:nyears,,,,1,simulation_number,drop=FALSE], c(1, 5), sum)
    om_bts_survey <- apply(survey_obs$rpws[1:nyears,,,,2,simulation_number,drop=FALSE], c(1, 5), sum)
    om_survey <- array(c(om_ll_survey, om_bts_survey), dim=c(nyears, 2), dimnames=list("time"=1:nyears, "fleet"=c("LLS", "BTS")))

    survey_df <- reshape2::melt(em_survey, value.name="em") %>%
        left_join(reshape2::melt(om_survey, value.name="om"), by=c("time", "fleet"))

    idx_plot <- ggplot(survey_df)+
        geom_line(aes(x=time, y=om), color="red", alpha=0.30)+
        # geom_point(aes(x=time, y=om), color="red")+
        geom_pointrange(aes(x=time, y=om, ymin=om-2.96*om/10, ymax=om+2.96*om/10), alpha=0.30, color="red", size=0.5)+
        geom_line(aes(x=time, y=em))+
        labs(x="Simulation Year", y="Index", title="Fits to Indices")+
        facet_wrap(~fleet, scales="free_y")+
        custom_theme

    # Plot catch fits
    em_catch <- aperm(report$PredCatch, c(2, 1, 3))[1:nyears,,1:2]
    dimnames(em_catch) <- list("time"=1:nyears, "fleet"=c("Fixed", "Trawl"))

    om_catch <- apply(catch[1:nyears,,,,,drop=FALSE], c(1, 5), sum)

    catch_df <- reshape2::melt(em_catch, value.name="em") %>%
        left_join(reshape2::melt(om_catch, value.name="om"), by=c("time", "fleet"))

    catch_plot <- ggplot(catch_df)+
        geom_line(aes(x=time, y=em))+
        geom_line(aes(x=time, y=om), color="red", alpha=0.30)+
        geom_point(aes(x=time, y=om), color="red")+
        labs(x="Simulation Year", y="Catch", title="Fits to Catch")+
        facet_wrap(~fleet)+
        custom_theme

    # Plot age comp fits
    em_srv_ac <- aperm(report$SrvIAA, c(2, 3, 4, 1, 5))
    em_srv_ac <- array(aperm(apply(em_srv_ac, c(1, 3, 5), \(x) x/sum(x)), c(2, 1, 3, 4)), dim=c(nyears, 30, 2, 1, 2))
    dimnames(em_srv_ac) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LLS", "BTS"))

    em_fish_ac <- aperm(report$CAA, c(2, 3, 4, 1, 5))
    em_fish_ac <- array(aperm(apply(em_fish_ac, c(1, 3, 5), \(x) x/sum(x)), c(2, 1, 3, 4)), dim=c(nyears, 30, 2, 1, 2))
    dimnames(em_fish_ac) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL"))

    em_ac <- array(c(em_fish_ac, em_srv_ac), dim=c(nyears, 30, 2, 1, 4))
    dimnames(em_ac) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL", "LLS", "BTS"))

    om_ac <- afscOM::subset_matrix(survey_obs$agg_acs, simulation_number, 6, drop=TRUE) 
    om_ac <- afscOM::subset_matrix(om_ac, 1:nyears, 1, drop=FALSE)
    om_ac <- array(aperm(apply(om_ac, c(1, 3, 5), \(x) x/sum(x)), c(2, 1, 3, 4)), dim=c(nyears, 30, 2, 1, 4))
    dimnames(om_ac) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL", "LLS", "BTS"))

    ac_df <- reshape2::melt(em_ac, value.name="em") %>% as_tibble() %>%
        left_join(reshape2::melt(om_ac, value.name="om"), by=c("time", "age", "sex", "region", "fleet"))

    ac_plot <- ggplot(ac_df)+
        geom_point(aes(x=om, y=em, color=age, shape=sex))+
        geom_abline(slope=1)+
        scale_y_continuous("Estimated", limits=c(0, 0.8))+
        scale_x_continuous("True", limits=c(0, 0.8))+
        scale_color_viridis_c()+
        labs(title="Fits to ACs")+
        facet_grid(fleet~sex)+
        custom_theme

    # Plot selectivity fits
    em_fish_sel <- aperm(report$fish_sel, c(2, 3, 4, 1, 5))
    dimnames(em_fish_sel) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL"))

    em_srv_sel <- aperm(report$srv_sel, c(2, 3, 4, 1, 5))
    dimnames(em_srv_sel) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LLS", "BTS"))

    em_sel <- array(c(em_fish_sel, em_srv_sel), dim=c(nyears, 30, 2, 1, 4))
    dimnames(em_sel) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL", "LLS", "BTS"))


    om_fish_sel <- afscOM::subset_matrix(dem_params$sel[1:nyears,,,,], 1, 4, drop=FALSE)
    dimnames(om_fish_sel) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL"))

    om_srv_sel <- afscOM::subset_matrix(dem_params$surv_sel[1:nyears,,,,], 1, 4, drop=FALSE)
    dimnames(om_srv_sel) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LLS", "BTS"))

    om_sel <- array(c(om_fish_sel, om_srv_sel), dim=c(nyears, 30, 2, 1, 4))
    dimnames(om_sel) <- list("time"=1:nyears, "age"=2:31, "sex"=c("F", "M"), "region"="Alaska", "fleet"=c("LL", "TWL", "LLS", "BTS"))

    sel_df <- reshape2::melt(em_sel, value.name="em") %>% as_tibble() %>%
        left_join(reshape2::melt(om_sel, value.name="om"), by=c("time", "age", "sex", "region", "fleet"))

    sel_plot <- ggplot(sel_df)+
        geom_point(aes(x=age, y=om, group=time), color="red", size=0.75)+
        geom_line(aes(x=age, y=om, group=time), color="red", alpha=0.30)+
        geom_line(aes(x=age, y=em, group=time))+
        labs(title="Fits to Selectivity")+
        facet_grid(fleet~sex)+
        custom_theme


    # Combine plots

    #  plot <- (ssb_plot+hr_plot+rec_plot)/(catch_plot+idx_plot)/(ac_plot)
    plot <- (ssb_plot+rec_plot+(hr_plot/idx_plot))/(ac_plot+sel_plot) + 
        plot_annotation(
            title = paste0(
                "OM: ", mse_obj$om$name, " | ", "HCR: ", mse_obj$mp$name
            ),
            subtitle = paste0(
                "Simulation Year: ", simulation_year, " | Simulation Seed: ", simulation_seed
            ),
            caption = paste0("Convergence Message:\n", convergence_message)
        )

    return(plot)
}

#' Plot Map of Alaska with NMFS Stat Areas Colored
#' 
#' Original code from Matt Cheng
plot_alaska_map <- function(regions = c("BS", "AI", "WGOA", "CGOA", "EGOA")){
    # Read in maps here
    west = rnaturalearth::ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
    west = sf::st_shift_longitude(west) # shift ongitude for plotting

    # Read in stat areas
    nmfs_areas = sf::read_sf(dsn = here::here("data", "NMFS_Stat_Areas"), layer = "Sablefish_Longline_Area")
    nmfs_areas = nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A")) %>% 
    mutate(         NAME = case_when(
        NAME == "Aleutian Islands" ~ "AI",
        NAME == "Bering Sea" ~ "BS",
        NAME == "Western Gulf of Alaska" ~ "WGOA",
        NAME == "Central Gulf of Alaska" ~ "CGOA",
        NAME == "West Yakutat" ~ "EGOA",
        NAME == "East Yakutat / Southeast Alaska" ~ "EGOA"
    ), NAME = factor(NAME, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>% 
    group_by(NAME) %>%
    summarise(geometry = st_union(geometry))

    # Coerce longline areas
    nmfs_areas = st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
    nmfs_areas = nmfs_areas %>% st_transform(4326) # transform to crs 4326
    nmfs_areas = st_shift_longitude(nmfs_areas) # shift longitude for plotting

    # Set colors
    true_colors = c("BS" = "#E31C39", "AI" = "#EA8115", "WGOA" = "#1C39E3", "CGOA" = "#30AF6C", "EGOA" = "#8115EA")
    alphas = rep(0, length(true_colors))
    names(alphas) = names(true_colors)
    alphas[names(alphas) %in% regions] = 0.55

    labels = data.frame(
        x = c(190, 180, 195, 206, 218),
        y = c(58, 50, 52, 55, 57.5),
        name = c("BS", "AI", "WGOA", "CGOA", "EGOA")
    )

    map <- ggplot() +
        geom_sf(data = nmfs_areas, aes(fill = NAME, alpha=NAME)) +
        geom_sf(data = west, lwd = 0.25, alpha = 1) +
        geom_text(data=labels, aes(x=x, y=y, label=name), size=6)+
        coord_sf(ylim = c(45.2, 70.5), xlim = c(165, 230)) + # Restrict Map Area
        scale_alpha_manual(values = alphas) +
        #   scale_fill_manual(values = colors) +
        guides(fill="none")+
        labs(x = "Longitude", y = "Latitude")+
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank(),
        )
    
    return(map)
}

set_hcr_colors <- function(hcrs){
    hcr_colors <- scales::hue_pal()(length(hcrs))
    hcr_colors[which(hcrs == "No Fishing")] <- "#000000"
    names(hcr_colors) <- hcrs
    return(hcr_colors)
}

set_hcr_colors2 <- function(hcrs){
    hcr_colors <- c(
        "F40" = "#E31C39",
        "F50" = "#EA8115",
        "B40/F50" = "#1C39E3",
        "No Fishing" = "#000000",
        "F40 +/- 5%" = "#30AF6C",
        "F40 +/- 10%" = "#8115EA",
        "15k Harvest Cap" = "#29C1D6",
        "25k Harvest Cap" = "#AC7A53",
        "Constant F50" = "#f3ac0c",
        "PFMC 40-10" = "#83c738",
        "British Columbia" = "#FA5CCC"
    )
    return(hcr_colors)
}


rank_colors <- c(
    "#D55E00",
    "#FF740A",
    "#FF8B33",
    "#FFA35C",
    "#FFBA85",
    "#AAAAAA",
    "#5CC3FF",
    "#33B4FF",
    "#0AA5FF",
    "#008EE0",
    "#0072B2"
)

custom_theme <- ggplot2::theme_bw()+ggplot2::theme(
    panel.spacing.y = ggplot2::unit(0.5, "cm"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size=14),
    axis.text = ggplot2::element_text(size=14),
    strip.text = ggplot2::element_text(size=14),
    legend.text = ggplot2::element_text(size=14),
    legend.position = "bottom"
)

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

plot_timeseries <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40", ylab=""){

    hcr1 <- as.character((data %>% pull(hcr) %>% unique)[1])

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- data %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- data %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(eval(rlang::parse_expr(v2))) %>% filter(time <= common)

    base_hcr_d <- data %>% filter(hcr == base_hcr)

    plot <- ggplot(data) + 
        geom_line(data = base_hcr_d, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        scale_fill_brewer(palette="Greys")+
        scale_color_manual(values=hcr_colors)+
        labs(x="Year", y=ylab)+
        coord_cartesian(expand=0)+
        guides(color=guide_legend(title="Management \n Strategy", nrow=2), fill="none")

    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(.data[[v3]]))+guides(fill="none")
    }
    
    return(plot+custom_theme)
}


