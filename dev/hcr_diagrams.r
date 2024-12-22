rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)
# library(afscOM)
# library(afscOM) # may work but not certain

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

devtools::load_all(afscOM_dir)

lapply(list.files("R", full.names = TRUE), source)

sable_om <- readRDS(file.path(here::here(), "data", "sablefish_om_big.RDS"))
assessment <- dget(file.path(here::here(), "data", "sablefish_assessment_2023.rdat"))
hist_recruits <- assessment$natage.female[,1]*2

dp_y <- afscOM::subset_dem_params(sable_om$dem_params, 64, d=1, drop=FALSE)
joint_selret <- calculate_joint_selret(
    sel = dp_y$sel,
    ret = dp_y$ret,
    prop_fs = c(0.80, 0.20)
)

spr40_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.40, 0.40)
)

spr50_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.50, 0.50)
)

spr4050_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.50, 0.40)
)

sprMSY_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.45, 0.45)
)

convert_to_catch <- function(ssb, f){
    std_sel <- joint_selret$sel[,,1,]/sum(joint_selret$sel[,,1,])
    catch <- afscOM::F_to_mu(f)*sum(ssb*std_sel)
    return(catch)
}

threshold_f_hcr <- function(ssb, FRPs, BRPs){
    return(threshold_f(ssb, f_min=FRPs[1], f_max=FRPs[2], lrp=BRPs[1], urp=BRPs[2]))
}

threshold_f_dep_hcr <- function(ssb, FRPs, BRPs, Bref){
    dep <- ssb/Bref
    return(threshold_f_hcr(dep, FRPs, BRPs))
}

threshold_cap_hcr <- function(ssb, FRPs, BRPs, cap){
    f1 <- threshold_f_hcr(ssb, FRPs, BRPs)
    std_sel <- joint_selret$sel[,,1,]/sum(joint_selret$sel[,,1,])

    catch <- afscOM::F_to_mu(f1)*sum(ssb*std_sel)
    if(catch > cap) catch <- cap

    mu <- catch/sum(ssb*std_sel)
    f2 <- afscOM::mu_to_F(mu)
    return(f2)
}



ssbs <- seq(0, 350, 1)

f40_hcr <- sapply(ssbs, \(x) threshold_f_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref)))
f50_hcr <- sapply(ssbs, \(x) threshold_f_hcr(x, FRPs=c(0, spr50_rp$Fref), BRPs=c(spr50_rp$Bref*0.05, spr50_rp$Bref)))
b40f50_hcr <- sapply(ssbs, \(x) threshold_f_hcr(x, FRPs=c(0, spr4050_rp$Fref), BRPs=c(spr4050_rp$Bref*0.05, spr4050_rp$Bref)))
pfmc4010_hcr <- sapply(ssbs, \(x) threshold_f_dep_hcr(x, FRPs=c(0, sprMSY_rp$Fref), BRPs=c(0.1, 0.40), Bref=sprMSY_rp$B0))
bcsable_hcr <- sapply(ssbs, \(x) threshold_f_dep_hcr(x, FRPs=c(0, 0.055), BRPs=c(0.40, 0.60), Bref=sprMSY_rp$Bref))
cap15_hcr <- sapply(ssbs, \(x) threshold_cap_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref), cap=15))
cap25_hcr <- sapply(ssbs, \(x) threshold_cap_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref), cap=25))
chrf50_hcr <- rep(spr50_rp$Fref, length(ssbs))

hcr_df <- data.frame(
        ssb=ssbs, 
        f40=f40_hcr,
        perc5=f40_hcr,
        perc10=f40_hcr, 
        f50=f50_hcr, 
        b40f50=b40f50_hcr, 
        pfmc4010=pfmc4010_hcr, 
        bcsable=bcsable_hcr, 
        chr50=chrf50_hcr,
        cap15=cap15_hcr,
        cap25=cap25_hcr
)

hcr_df_long <- hcr_df %>% pivot_longer(2:ncol(hcr_df), names_to="hcr", values_to="F") %>%
    mutate(
        hcr=factor(
            hcr, 
            levels=c("f40", "f50", "b40f50", "perc5", "perc10", "cap15", "cap25", "chr50", "pfmc4010", "bcsable"),
            labels=c("F40", "F50", "B40/F50", "F40 +/- 5%", "F40- +/- 10%", "15k Harvest Cap", "25k Harvest Cap", "Constant F50", "PFMC 40-10", "British Columbia")
        )
    ) %>%
    rowwise() %>%
    mutate(
        catch = convert_to_catch(ssb, F)
    )


ggplot(hcr_df_long)+
    geom_line(aes(x=ssb, y=F, color=hcr, linetype="F"), linewidth=1)+
    geom_line(aes(x=ssb, y=catch/267, color=hcr, linetype="Catch"), linewidth=1)+
    geom_vline(xintercept = spr40_rp$Bref, linetype="dashed", color="grey")+
    geom_vline(xintercept = spr50_rp$B0, linetype="dashed", color="grey")+
    annotate("text", x=spr40_rp$Bref-40, y=0.11, label="B[40]", parse=TRUE, size=6)+
    annotate("text", x=spr40_rp$B0-40, y=0.11, label="B[100]", parse=TRUE, size=6)+
    scale_x_continuous("Spawning Stock Biomass (SSB)")+
    scale_y_continuous("F", limits=c(0, 0.12), sec.axis = sec_axis(transform = ~.*267, name="Catch (1k mt)"))+
    scale_linetype_manual(values=c("F"="solid", "Catch"="dotdash"))+
    coord_cartesian(expand=0)+
    facet_wrap(~hcr, nrow=2)+
    guides(color="none")+
    custom_theme+
    theme(panel.grid = element_blank())
ggsave(file.path(here::here(), "figures", "hcr_diagrams.jpeg"), width=15, units="in")

catch_plot <- ggplot(hcr_df_long)+
    geom_line(aes(x=ssb, y=catch, color=hcr), linewidth=1)+
    geom_vline(xintercept = spr40_rp$Bref, linetype="dashed", color="grey")+
    geom_vline(xintercept = spr50_rp$B0, linetype="dashed", color="grey")+
    annotate("text", x=spr40_rp$Bref-40, y=28, label="B[40]", parse=TRUE, size=6)+
    annotate("text", x=spr40_rp$B0-40, y=28, label="B[0]", parse=TRUE, size=6)+
    scale_y_continuous("Catch (1000s mt)", limits=c(0, 32))+
    coord_cartesian(expand=0)+
    facet_wrap(~hcr, nrow=2)+
    custom_theme
ggsave(file.path(here::here(), "figures", "hcr_catch_diagrams.jpeg"), width=15, units="in")

f_plot <- ggplot(hcr_df_long)+
    geom_line(aes(x=ssb, y=F, color=hcr), linewidth=1)+
    geom_vline(xintercept = spr40_rp$Bref, linetype="dashed", color="grey")+
    geom_vline(xintercept = spr50_rp$B0, linetype="dashed", color="grey")+
    annotate("text", x=spr40_rp$Bref-40, y=0.11, label="B[40]", parse=TRUE, size=6)+
    annotate("text", x=spr40_rp$B0-40, y=0.11, label="B[0]", parse=TRUE, size=6)+
    scale_y_continuous("Fishing Mortality", limits=c(0, 0.12))+
    coord_cartesian(expand=0)+
    facet_wrap(~hcr, nrow=2)+
    custom_theme
ggsave(file.path(here::here(), "figures", "hcr_f_diagrams.jpeg"), width=15, units="in")

library(patchwork)

f_plot / catch_plot + plot_layout(guides="collect") & theme(legend.position = "none")
