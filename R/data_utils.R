reformat_ggdist_long <- function(data){

    data_long <- data %>% pivot_longer(-c(time, .width, .point, .interval))
    medians <- data_long %>% filter(!(grepl("lower",name) | grepl("upper",name)))
    lowers <- data_long %>% filter(grepl("lower", name)) %>% mutate(name=unlist(lapply(str_split(name, fixed(".")), `[[`, 1)))
    uppers <- data_long %>% filter(grepl("upper", name)) %>% mutate(name=unlist(lapply(str_split(name, fixed(".")), `[[`, 1)))

    return(
        medians %>% 
            left_join(lowers, by=c("time", ".width", ".point", ".interval", "name")) %>% 
            left_join(uppers, by=c("time", ".width", ".point", ".interval", "name")) %>%
            rename("median"="value.x", "lower"="value.y", "upper"="value")
    )
    
}