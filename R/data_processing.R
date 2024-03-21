create_summary_biomass_df <- function(naa, waa, mat){
    reshape2::melt(naa, value.name="naa") %>% 
        as_tibble() %>%
        #filter(sex=="F") %>%
        dplyr::left_join(
            reshape2::melt(waa, value.name="weight"), 
            by=c("time", "age", "sex")
        ) %>%
        dplyr::left_join(
            reshape2::melt(mat, value.name="maturity"), 
            by=c("time", "age", "sex")
        ) %>%
        drop_na() %>%
        mutate(
            biomass = naa*weight,
            spbio = naa*weight*maturity
        ) %>%
        select(-c(starts_with("region"), weight, maturity, naa))
}

create_summary_catch_df <- function(caa, out_f){
    reshape2::melt(caa, value.name="caa") %>% as_tibble() %>%
        left_join(
            reshape2::melt(out_f, value.name = "F"), by=c("time", "sims")
        ) %>%
        select(-c(starts_with("region"), Var2, Var3)) %>%
        group_by(time, age, sex, sims) %>%
        summarise(
            catch = sum(caa),
            f=mean(F)
        ) %>%
        ungroup()
}

create_biomass_catch_summary_df <- function(bio_summ, catch_summ){
    bio_summ %>%
        left_join(
            catch_summ,
            by=c("time", "age", "sex", "sim"="sims")
        ) %>%
        group_by(time, sim) %>%
        summarise(
            tot_bio = sum(biomass),
            tot_spbio = sum(spbio[sex == "F"]),
            catch=sum(catch),
            F=f
        ) %>%
        group_by(time) %>%
        median_qi(tot_bio, tot_spbio, catch, F, .width=c(0.5, 0.95), .simple_names=FALSE) %>%
        reformat_ggdist_long %>%
        mutate(
            name = factor(name, levels=c("tot_spbio", "catch", "tot_bio", "F"), labels=c("Spawning Biomass (mt)", "Catch (mt)", "Total Biomass (mt)", "Summary F"))
        )
}
