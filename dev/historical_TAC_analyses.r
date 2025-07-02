library(tidyverse)

historical_data <- read_csv(
    file.path("data", "sablefish_historical_regional_TACs.csv")
)

historical_data_long <- historical_data %>% slice(c(2, 3, 6, 8, 9)) %>%
    setNames(c("Region", paste0(rep(c(2024:1986), each=3), "_", rep(c("OFL", "ABC", "TAC"), 3)))) %>%
    pivot_longer(
        2:ncol(.), 
        names_to="quantity",
        values_to="value"
    ) %>%
    separate(
        quantity, 
        into=c("Year", "quantity"), 
        sep="_"
    ) %>%
    mutate(
        value = as.numeric(gsub(",", "", value)),
        Year = as.integer(Year),
        quantity = factor(quantity, levels=c("OFL", "ABC", "TAC"), ordered=TRUE),
        Region = factor(Region, levels=c("BS", "AI", "WGOA", "CGOA", 'EGOA'))
    ) %>%
    arrange(Year, Region, quantity)

hist_fleet_tacs <- historical_data_long %>% filter(Year > 1999, quantity != "OFL") %>%
    pivot_wider(
        names_from=quantity,
        values_from=value
    ) %>%
    mutate(
        reduction = TAC/ABC,
        fixed_split = case_when(
            Region == "BS" ~ 0.5,
            Region == "AI" ~ 0.75,
            Region == "WGOA" ~ 0.8,
            Region == "CGOA" ~ 0.8,
            Region == "EGOA" ~ 0.95
        ),
        trawl_split = case_when(
            Region == "BS" ~ 0.5,
            Region == "AI" ~ 0.25,
            Region == "WGOA" ~ 0.2,
            Region == "CGOA" ~ 0.2,
            Region == "EGOA" ~ 0.05
        ),
        fixed_tac = TAC* fixed_split,
        trawl_tac = TAC * trawl_split
    ) %>%
    select(-fixed_split, -trawl_split)

hist_fleet_tacs %>%
    ggplot(aes(x=Year, y=reduction, color=Region)) +
        geom_line() +
        geom_point() +
        labs(
            title="Historical Sablefish ABC-TAC Reduction by Region",
            x="Year",
            y="TAC (mt)",
            color="Region"
        ) +
        scale_y_continuous(labels=scales::percent) +
        coord_cartesian(ylim=c(0, 1.2)) +
        custom_theme
ggsave(file.path(here::here(), "figures", "historical_management", "sablefish_abc_tac_reduction_by_region.jpeg"))

hist_fleet_tacs %>% select(Region, Year, fixed_tac, trawl_tac) %>%
    rename(Fixed = fixed_tac, Trawl = trawl_tac) %>%
    pivot_longer(
        cols=c(Fixed, Trawl),
        names_to="Fleet",
        values_to="TAC"
    ) %>%
    mutate(Fleet = factor(Fleet, labels=c("Fixed Gear", "Trawl"))) %>%

    ggplot(aes(x=Year, y=TAC, color=Region, linetype=Fleet)) +
        geom_line() +
        geom_point() +
        labs(
            title="Historical Sablefish TACs by Region",
            x="Year",
            y="TAC (mt)",
            color="Region"
        ) +
        custom_theme


assessment <- dget(file.path(here::here(), "data/spatial_sablefsh_inputs.rdat"))

catch_data <- reshape2::melt(assessment$catch) %>%
    setNames(c("Region", "Year", "Fleet", "Landings")) %>%
    mutate(
        Region = factor(Region, labels=c("BS", "AI", "WGOA", "CGOA", 'EGOA')),
        Year = rep(rep(assessment$years, each=5), 2),
        Fleet = factor(Fleet, labels=c("Fixed Gear", "Trawl")),
        Landings = Landings*1000
    )

regflt_tac_landings <- hist_fleet_tacs %>% select(Region, Year, fixed_tac, trawl_tac) %>%
    rename(Fixed = fixed_tac, Trawl = trawl_tac) %>%
    pivot_longer(
        cols=c(Fixed, Trawl),
        names_to="Fleet",
        values_to="TAC"
    ) %>%
    mutate(Fleet = factor(Fleet, labels=c("Fixed Gear", "Trawl"))) %>%
    filter(Year < 2022) %>%
    left_join(
        catch_data %>% filter(Year > 1999)
    ) %>%
    mutate(
        utilization = Landings / TAC,
    )

ggplot(regflt_tac_landings, aes(x=Year, y=utilization, color=Region, linetype=Fleet)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept=1, linetype="dashed", color="black") +
    labs(
        title="Historical Sablefish TAC Utilization by Region and Fleet",
        x="Year",
        y="Utilization (Landings / TAC)",
        color="Region"
    ) +
    facet_wrap(~Region)+
    custom_theme+
    scale_y_continuous(labels=scales::percent)+
    coord_cartesian(ylim=c(0, 2.0)) 
ggsave(file.path(here::here(), "figures", "historical_management", "sablefish_tac_utilization_by_region_fleet.jpeg"))

regflt_tac_landings %>% group_by(Year, Region) %>%
    summarise(
        TAC = sum(TAC),
        Landings = sum(Landings),
        utilization = Landings/TAC
    ) %>%

    ggplot(aes(x=Year, y=utilization, color=Region)) +
        geom_line() +
        geom_point() +
        geom_hline(yintercept=1, linetype="dashed", color="black") +
        labs(
            title="Historical Sablefish TAC Utilization by Region",
            x="Year",
            y="Utilization (Landings / TAC)",
            color="Region"
        ) +
        custom_theme+
        scale_y_continuous(labels=scales::percent) +
        coord_cartesian(ylim=c(0, 1.5))
ggsave(file.path(here::here(), "figures", "historical_management", "sablefish_tac_utilization_by_region.jpeg"))


regflt_tac_landings %>% group_by(Year) %>%
    summarise(
        TAC = sum(TAC),
        Landings = sum(Landings),
        utilization = Landings/TAC
    ) %>%

    ggplot(aes(x=Year, y=utilization)) +
        geom_line() +
        geom_point() +
        geom_hline(yintercept=1, linetype="dashed", color="black") +
        labs(
            title="Historical Sablefish TAC Utilization",
            x="Year",
            y="Utilization (Landings / TAC)"
        ) +
        custom_theme+
        scale_y_continuous(labels=scales::percent) +
        coord_cartesian(ylim=c(0, 1.5))
ggsave(file.path(here::here(), "figures", "historical_management", "sablefish_tac_utilization_alaska.jpeg"))





regflt_tac_landings %>%

    ggplot(aes(x=TAC, y=utilization, color=Region, shape=Fleet))+
        geom_smooth(method="lm", se=TRUE, linetype="dashed", color="black")+
        geom_point(size=3) +
        geom_hline(yintercept=1, linetype="dashed", color="black") +
        coord_cartesian(ylim=c(0, 2.0))+
        facet_grid(Region ~ Fleet, scales="free") +
        labs(
            title="Historical Sablefish TAC vs Utilization by Region and Fleet",
            x="TAC (mt)",
            y="Utilization (Landings / TAC)",
            color="Region",
            shape="Fleet"
        )+
        custom_theme
ggsave(file.path(here::here(), "figures", "historical_management", "sablefish_tac_vs_utilization_by_region_fleet.jpeg"))
        
regflt_tac_landings %>%

    ggplot(aes(x=TAC, y=Landings, color=Region, shape=Fleet))+
        geom_smooth(method="lm", se=TRUE, linetype="dashed", color="black")+
        geom_point(size=3) +
        # geom_hline(yintercept=1, linetype="dashed", color="black") +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
        # coord_cartesian(ylim=c(0, 2.0))+
        facet_grid(Region ~ Fleet, scales="free") +
        labs(
            title="Historical Sablefish TAC vs Landings by Region and Fleet",
            x="TAC (mt)",
            y="Landings (mt)",
            color="Region",
            shape="Fleet"
        )+
        custom_theme
ggsave(file.path(here::here(), "figures", "historical_management", "sablefish_tac_vs_landings_by_region_fleet.jpeg"))



regflt_tac_landings %>% group_by(Region, Fleet) %>%
    summarise(
        utilization = median(utilization)
    ) %>%
    pivot_wider(
        names_from=Fleet,
        values_from=utilization
    )

hist_fleet_tacs %>% select(Region, Year, ABC, TAC, reduction) %>%
    group_by(Region) %>%
    summarise(
        reduction = min(reduction)
    ) %>%
    pivot_wider(
        names_from=Region,
        values_from=reduction
    )

library(betareg)

bs_trawl_util <- regflt_tac_landings %>%
    filter(Region == "BS", Fleet == "Trawl") %>%
    select(TAC, utilization) %>%
    mutate(
        utilization = ifelse(utilization > 1, 1, utilization)  # Cap utilization at 1
    )

reshape2::melt(sable_om$dem_params$movement[,,200,,1]) %>%
    setNames(c("From", "To", "Age", "Movement")) %>%
    filter(Age %in% c(2, 13, 25)) %>%
    ggplot(aes(x=From, y=To, fill=Movement)) +
        geom_tile() +
        scale_fill_gradient(low="white", high="red", limits=c(0, 1)) +
        facet_wrap(~Age, nrow=3) +
        labs(
            title="Movement Matrix",
            x="From Region",
            y="To Region",
            fill="Movement"
        ) +
        custom_theme
ggsave("~/Desktop/sablefish_movement_matrix.jpeg", height=9, width=5, units="in")

reshape2::melt(hist_recruits) %>%
    setNames(c("Year", "Region", "Recruits")) %>%
    mutate(
        Region = factor(Region, labels=c("BS", "AI", "WGOA", "CGOA", 'EGOA')),
        Year = Year+1960
    ) %>%
    ggplot(aes(x=Year, y=Recruits, color=Region)) +
        geom_line() +
        geom_point() +
        labs(
            title="Historical Sablefish Recruits by Region",
            x="Year",
            y="Recruits (thousands)",
            color="Region"
        ) +
        facet_wrap(~Region) +
        custom_theme
ggsave(file.path(here::here(), "figures", "sablefish_recruits_by_region.jpeg"))

reshape2::melt(hist_recruits) %>%
    setNames(c("Year", "Region", "Recruits")) %>%
    mutate(
        Region = factor(Region, labels=c("BS", "AI", "WGOA", "CGOA", 'EGOA')),
        Year = Year+1960
    ) %>%
    group_by(Year) %>%
    mutate(total_recruits = sum(Recruits)) %>%
    ungroup() %>%
    mutate(
        proportion = Recruits / total_recruits
    ) %>%
    ggplot(aes(x=Year, y=proportion, color=Region)) +
        geom_line() +
        geom_point() +
        labs(
            title="Historical Proportion of Sablefish Recruits by Region",
            x="Year",
            y="Proportion of Total Recruits",
            color="Region"
        ) +
        facet_wrap(~Region) +
        custom_theme
ggsave(file.path(here::here(), "figures", "proportion_recruits_by_region.jpeg"))

cov(hist_recruits)

apply(t(apply(hist_recruits, 1, \(x) x/sum(x))), 2, mean)

cor(hist_recruits)


historical_data %>% slice(c(2, 3, 6, 8, 9, 10)) %>%
    setNames(c("Region", paste0(rep(c(2024:1986), each=3), "_", rep(c("OFL", "ABC", "TAC"), 3)))) %>%
    mutate(Region=ifelse(is.na(Region), "AK", Region)) %>%
    pivot_longer(
        2:ncol(.), 
        names_to="quantity",
        values_to="value"
    ) %>%
    separate(
        quantity, 
        into=c("Year", "quantity"), 
        sep="_"
    ) %>%
    mutate(
        value = as.numeric(gsub(",", "", value)),
        Year = as.integer(Year),
        quantity = factor(quantity, levels=c("OFL", "ABC", "TAC"), ordered=TRUE),
        Region = factor(Region, levels=c("BS", "AI", "WGOA", "CGOA", 'EGOA', "AK"))
    ) %>%
    arrange(Year, Region, quantity) %>%
    filter(Region == "AK", quantity == "ABC", Year > 2014) %>%
    pull(value) %>% summary
