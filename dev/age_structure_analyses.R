#' 9. Compute age structure summary statistics
#' These statistics include Shannon diversity (Shannon and Weaver 1949),
#' average age (van Deurs et al. 2023), proportion mature (van Deurs et al. 2023),
#' proportion fully mature (Goethel, pers. comm), and ABI (Griffiths et al. 2023).

get_om_data <- function(fname){
    data <- readRDS(fname)
    naa <- data$naa
    dem_params <- data$dem_params

    total_naa <- apply(naa, c(1, 2, 5), sum)
    total_naa <- total_naa[1:1000,,]

    female_naa <- naa[1:1000,,1,1,]
    female_baa <- array(apply(female_naa, 3, \(s) s*dem_params$waa[,,1,1]), dim=c(nyears, 30, 100))
    female_sbaa <- array(apply(female_naa, 3, \(s) s*dem_params$waa[,,1,1]*dem_params$mat[,,1,1]), dim=c(nyears, 30, 100))
    female_mat <- dem_params$mat[1:1000,,1,1]

    male_naa <- naa[1:1000,,2,1,]
    male_baa <- array(apply(male_naa, 3, \(s) s*dem_params$waa[,,2,1]), dim=c(nyears, 30, 100))

    total_baa <- female_baa+male_baa

    recruits  <- apply(total_naa[,1,], 1, mean)
    return(afscOM::listN(female_naa, total_baa, female_sbaa, female_mat, recruits))
}


generate_indices <- function(naa, mat){
    avg_agestruct <- apply(naa, c(1, 2), mean)
    avg_agestruct_equil <- apply(avg_agestruct[500:1000,], 2, mean)

    shan_div  <- apply(apply(naa, c(1, 3), shannon_diversity), 1, mean)
    amean     <- apply(apply(naa, c(1, 3), average_age, ages=2:31), 1, mean)
    pmat      <- apply(sapply(1:100, \(s) sapply(1:nyears, \(x) prop_mature(naa=naa[x,,s], mat=mat[x,]))), 1, mean)
    pmat_full <- apply(sapply(1:100, \(s) sapply(1:nyears, \(x) prop_fully_mature(naa=naa[x,,s], mat=mat[x,]))), 1, mean)
    abi0      <- apply(apply(naa, c(1, 3), abi, ref_naa=avg_agestruct_equil), 1, mean)
    

    as_df <- data.frame(year=1:nyears, shannon=shan_div, avgage=amean, propmat=pmat, propmatf=pmat_full, abi=abi0) %>%
        pivot_longer(-c(year), names_to="metric", values_to="value")

    return(as_df)
}

om_00 <- get_om_data("data/om_f00.RDS")

#' Tests with numbers at age distribution
as_df_naa <- generate_indices(om_00$female_naa, om_00$female_mat)

ggplot(as_df_naa %>% filter(metric != "rec"))+
    geom_line(aes(x=year, y=value, color=metric))+
    facet_wrap(~metric, scales="free_y")

#' Tests with total biomass at age distribution
as_df_baa <- generate_indices(om_00$total_baa, om_00$female_mat)

ggplot(as_df_baa %>% filter(metric != "rec"))+
    geom_line(aes(x=year, y=value, color=metric))+
    facet_wrap(~metric, scales="free_y")

#' Tests with spawning biomass at age distribution
as_df_sbaa <- generate_indices(om_00$female_sbaa, om_00$female_mat)

ggplot(as_df_sbaa %>% filter(metric != "rec"))+
    geom_line(aes(x=year, y=value, color=metric))+
    facet_wrap(~metric, scales="free_y")

#' Combine plot

library(ggh4x)
age_struct_summ_00 <- generate_age_struct_summ("data/om_f00.RDS")

pos_scales <- list(
    scale_y_continuous(limits=c(0, 4.1), breaks=seq(0, 4, 1), sec.axis = sec_axis(~ . *25, breaks=seq(0, 100, 25))),
    scale_y_continuous(limits=c(0, 20),  breaks=seq(0, 20, 5), sec.axis = sec_axis(~ . *5, breaks=seq(0, 100, 25))),
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.25), sec.axis = sec_axis(~ . *100, breaks=seq(0, 100, 25), name = "Recruitment (millions)")),
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.25), sec.axis = sec_axis(~ . *100, breaks=seq(0, 100, 25))),
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.25), sec.axis = sec_axis(~ . *100, breaks=seq(0, 100, 25)))
)

ggplot(age_struct_summ_00 %>% filter(year < 64))+
    geom_line(aes(x=year, y=value, color=metric, linetype=unit), size=0.85)+
    geom_line(aes(x=year, y=rec), color="black", alpha=0.5)+
    geom_vline(xintercept=20, color="grey")+
    geom_vline(xintercept=57, color="grey")+
    scale_x_continuous(limits=c(1, 65), breaks=seq(10, 64, 20), labels=seq(1970, 2023, 20), expand=c(0, 0))+
    #facet_grid(rows=vars(metric), cols=vars(unit), scales="free_y")+
    facet_wrap(~metric, scales="free_y")+
    facetted_pos_scales(y=pos_scales)+
    guides(color="none")+
    labs(x="Year", y="Value")+
    theme_bw()+
    theme(
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=24),
        strip.text = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)
    )
ggsave("~/Desktop/unit_comparison2.jpg")

generate_age_struct_summ <- function(om_fname){
    om <- get_om_data(om_fname)
    as_df_naa <- generate_indices(om$female_naa, om$female_mat)
    as_df_baa <- generate_indices(om$total_baa, om$female_mat)
    as_df_sbaa <- generate_indices(om$female_sbaa, om$female_mat)

    recruits <- om$recruits
    rec_df <- data.frame(year=rep(rep(1:nyears, each=5), 3), rec=c(rbind(recruits/25, recruits/5, recruits/100, recruits/100, recruits/100)))

    age_struct_summ <- bind_rows(
        as_df_naa %>% mutate(unit="numbers"),
        as_df_baa %>% mutate(unit="biomass"),
        as_df_sbaa %>% mutate(unit="spawning")
    ) %>% 
        mutate(
            rec=rec_df$rec,
            metric = factor(metric, levels=c("shannon", "avgage", "propmat", "propmatf", "abi"), labels=c("Diversity", "Average Age", "Prop Mature", "Fully Mature", "ABI0")),
            unit   = factor(unit, levels=c("numbers", "biomass", "spawning"), labels=c("Numbers-at-age", "Biomass-", "Spawning Biomass-"))
        )

    return(age_struct_summ)
}

compute_reference_levels <- function(age_struct_summ){
    age_struct_summ %>% 
        filter(year >= 500) %>% 
        group_by(metric, unit) %>% 
        summarise(
            equil_mean=mean(value)
        )
}

ref_levels_00 <- compute_reference_levels(age_struct_summ)
standardize_to_reference <- function(age_struct_summ, ref_levels){
    age_struct_summ %>%
        left_join(ref_levels, by=c("metric", "unit")) %>%
        mutate(value = value/equil_mean) %>%
        select(-c(equil_mean))
}

rel_age_struct_summ <- standardize_to_reference(age_struct_summ, ref_levels_00)
rel_age_struct_summ$rec <- rep(rep(recruits/100, each=5), 3)

ggplot(rel_age_struct_summ %>% filter(year < 64))+
    geom_line(aes(x=year, y=value, color=metric), size=0.85)+
    geom_line(aes(x=year, y=rec), color="black", alpha=0.5)+
    geom_vline(xintercept=20, color="grey")+
    geom_vline(xintercept=57, color="grey")+
    scale_x_continuous(limits=c(1, 65), breaks=seq(10, 64, 20), labels=seq(1970, 2023, 20), expand=c(0, 0))+
    facet_grid(rows=vars(metric), cols=vars(unit), scales="free_y")+
    #facetted_pos_scales(y=pos_scales)+
    guides(color="none")+
    labs(x="Year", y="Value")+
    theme_bw()+
    theme(
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=24),
        strip.text = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

as_summ_00 <- generate_age_struct_summ("data/om_f00.RDS")
as_summ_35 <- generate_age_struct_summ("data/om_f35.RDS")
as_summ_40 <- generate_age_struct_summ("data/om_f40.RDS")

age_struct_summ <- bind_rows(
    as_summ_00 %>% mutate(f="F0"),
    as_summ_35 %>% mutate(f="F35"),
    as_summ_35 %>% mutate(f="F40")
)

pos_scales <- list(
    scale_y_continuous(limits=c(0, 4.1), sec.axis = sec_axis(~ . *25, breaks=seq(0, 100, 25))),
    scale_y_continuous(limits=c(0, 20),  sec.axis = sec_axis(~ . *5, breaks=seq(0, 100, 25))),
    scale_y_continuous(limits=c(0, 1.1), sec.axis = sec_axis(~ . *100, breaks=seq(0, 100, 25), name = "Recruitment (millions)")),
    scale_y_continuous(limits=c(0, 1.1), sec.axis = sec_axis(~ . *100, breaks=seq(0, 100, 25))),
    scale_y_continuous(limits=c(0, 2.0), sec.axis = sec_axis(~ . *50, breaks=seq(0, 100, 25)))
)

ggplot(age_struct_summ %>% filter(year < 64))+
    geom_line(aes(x=year, y=value, color=metric, group=unit, linetype=unit), size=0.85)+
    geom_line(aes(x=year, y=rec), color="black", alpha=0.5)+
    geom_vline(xintercept=20, color="grey")+
    geom_vline(xintercept=57, color="grey")+
    scale_x_continuous(limits=c(1, 65), breaks=seq(10, 64, 20), labels=seq(1970, 2023, 20), expand=c(0, 0))+
    facet_grid(rows=vars(metric), cols=vars(f), scales="free_y")+
    facetted_pos_scales(y=pos_scales)+
    guides(color="none")+
    labs(x="Year", y="Value")+
    theme_bw()+
    theme(
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=24),
        strip.text = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

f00rl <- compute_reference_levels(as_summ_00)
f35rl <- compute_reference_levels(as_summ_35)
f40rl <- compute_reference_levels(as_summ_40)

age_struct_summ <- bind_rows(
    standardize_to_reference(as_summ_00, f00rl) %>% mutate(f="F0"),
    standardize_to_reference(as_summ_35, f35rl) %>% mutate(f="F35"),
    standardize_to_reference(as_summ_40, f40rl) %>% mutate(f="F40")
) %>%
    mutate(
        rec = rep(rep(rep(recruits/100, each=5), 3), 3)
    )

ggplot(age_struct_summ %>% filter(year < 64))+
    geom_line(aes(x=year, y=value, color=metric, group=unit, linetype=unit), size=0.85)+
    geom_line(aes(x=year, y=rec), color="black", alpha=0.5)+
    geom_vline(xintercept=20, color="grey")+
    geom_vline(xintercept=57, color="grey")+
    geom_hline(yintercept=1, color="grey")+
    scale_y_continuous(limits=c(0, 1.5), breaks=seq(0, 1.5, 0.5))+
    scale_x_continuous(limits=c(1, 65), breaks=seq(10, 64, 20), labels=seq(1970, 2023, 20), expand=c(0, 0))+
    facet_grid(rows=vars(metric), cols=vars(f), scales="free_y")+
    #facetted_pos_scales(y=pos_scales)+
    guides(color="none")+
    labs(x="Year", y="Value")+
    theme_bw()+
    theme(
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.title = element_text(size=24),
        strip.text = element_text(size=12),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing.y = unit(0.5, "cm")
    )
ggsave("~/Desktop/reference_levels.jpg")

#######

metrics <- as_summ_00 %>% pull(metric) %>% unique %>% levels
units   <- as_summ_00 %>% pull(unit) %>% unique %>% levels

m <- metrics[1]
u <- units[3]

recs <- as_summ_00 %>% filter(metric == m & unit == u & year < 64) %>% pull(rec)
vals <- as_summ_00 %>% filter(metric == m & unit == u & year < 64) %>% pull(value)

pdf(file="~/Desktop/ccfs.pdf")
par(mfrow=c(3, 3))
for(m in metrics[c(1,2,5)]){
    for(u in units){
        recs <- as_summ_00 %>% filter(metric == m & unit == u & year < 64) %>% pull(rec)
        vals <- as_summ_00 %>% filter(metric == m & unit == u & year < 64) %>% pull(value)
        title <- paste0(u, " ", m)
        ccf(recs, vals, main=title)
    }
}
dev.off()

####

om_f00 <- get_om_data("data/om_f00.RDS")
om_f35 <- get_om_data("data/om_f35.RDS")
om_f40 <- get_om_data("data/om_f40.RDS")

par(mfrow=c(2, 2))
barplot(apply(om_f00$female_naa[1000,,], 1, mean)/sum(apply(om_f00$female_naa[1000,,], 1, mean)), ylim=c(0, 0.15), main=expression(F["0"]))
barplot(apply(om_f35$female_naa[1000,,], 1, mean)/sum(apply(om_f35$female_naa[1000,,], 1, mean)), ylim=c(0, 0.15), main=expression(F["35"]))
barplot(apply(om_f40$female_naa[1000,,], 1, mean)/sum(apply(om_f40$female_naa[1000,,], 1, mean)), ylim=c(0, 0.15), main=expression(F["40"]))






# Shannon diversity Plot
ggplot(as_df %>% filter(year < 65 & metric %in% c("rec", "shannon")) %>% pivot_wider(names_from="metric", values_from="value"))+
    geom_line(aes(x=year, y=shannon), color="red", size=1)+
    geom_line(aes(x=year, y=rec/25), color="black", size=1)+
    scale_x_continuous(breaks=seq(1, 64, 20), labels=seq(1960, 2020, 20), expand=c(0, 0))+
    scale_y_continuous(limits=c(0, 4), breaks=seq(0, 4, 0.5), sec.axis=sec_axis(~ . * 25, breaks=seq(0, 100, 25), name="Recruitment (millions)"))+
    labs(x="Year", y="Shannon Diversity (H')", title="Shannon Diversity", subtitle="Shannon and Weaver, 1949")+
    theme_bw()+
    theme(
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=24)
    )
ggsave("~/Desktop/shannon.jpg")

# Proportion mature Plot
ggplot(as_df %>% filter(year < 65 & metric %in% c("rec", "avgage")) %>% pivot_wider(names_from="metric", values_from="value"))+
    geom_line(aes(x=year, y=avgage), color="red", size=1)+
    geom_line(aes(x=year, y=rec/8.333), color="black", size=1)+
    scale_x_continuous(breaks=seq(1, 64, 20), labels=seq(1960, 2020, 20), expand=c(0, 0))+
    scale_y_continuous(limits=c(0, 12), breaks=seq(0, 12, 2), sec.axis=sec_axis(~ . * 8.333, breaks=seq(0, 100, 25), name="Recruitment (millions)"))+
    labs(x="Year", y="Average Age", title="Average Age")+
    theme_bw()+
    theme(
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=24)
    )
ggsave("~/Desktop/avgage.jpg")

# Proportion mature Plot
ggplot(as_df %>% filter(year < 65 & metric %in% c("rec", "propmat")) %>% pivot_wider(names_from="metric", values_from="value"))+
    geom_line(aes(x=year, y=propmat), color="red", size=1)+
    geom_line(aes(x=year, y=rec/100), color="black", size=1)+
    scale_x_continuous(breaks=seq(1, 64, 20), labels=seq(1960, 2020, 20), expand=c(0, 0))+
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.25), sec.axis=sec_axis(~ . * 100, breaks=seq(0, 100, 25), name="Recruitment (millions)"))+
    labs(x="Year", y="Proportion Mature", title="Proportion Mature")+
    theme_bw()+
    theme(
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=24)
    )
ggsave("~/Desktop/propmat.jpg")

# Proportion fully mature Plot
ggplot(as_df %>% filter(year < 65 & metric %in% c("rec", "propmatf")) %>% pivot_wider(names_from="metric", values_from="value"))+
    geom_line(aes(x=year, y=propmatf), color="red", size=1)+
    geom_line(aes(x=year, y=rec/100), color="black", size=1)+
    scale_x_continuous(breaks=seq(1, 64, 20), labels=seq(1960, 2020, 20), expand=c(0, 0))+
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.25), sec.axis=sec_axis(~ . * 100, breaks=seq(0, 100, 25), name="Recruitment (millions)"))+
    labs(x="Year", y="Proportion Fully Mature", title="Proportion Fully Mature")+
    theme_bw()+
    theme(
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=24)
    )
ggsave("~/Desktop/propmat_full.jpg")

#' ABI Plot
ggplot(as_df %>% filter(year < 65 & metric %in% c("abi", "rec")) %>% pivot_wider(names_from="metric", values_from="value"))+
    geom_line(aes(x=year, y=abi), color="red", size=1)+
    geom_line(aes(x=year, y=rec/100), color="black", size=1)+
    geom_hline(yintercept=1, linetype="longdash", color="red")+
    scale_x_continuous(breaks=seq(1, 64, 20), labels=seq(1960, 2020, 20), expand=c(0, 0))+
    scale_y_continuous(limits=c(0, 1.1), breaks=seq(0, 1, 0.25), sec.axis=sec_axis(~ . * 100, breaks=seq(0, 100, 25), name="Recruitment (millions)"))+
    labs(x="Year", y=expression(ABI["0"]), title="ABI", subtitle="Griffiths et al. 2023")+
    theme_bw()+
    theme(
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        plot.title = element_text(size=24)
    )
ggsave("~/Desktop/abi.jpg")    


as_df %>% 
    filter(year < 65 & metric %in% c("abi", "rec")) %>% 
    pivot_wider(names_from="metric", values_from="value")
