OFL <- c(34100, 33200, 35900, 25700, 22800, 45600, 27800, 24700, 21500, 20700, 26100, 28900, 30800, 25400, 25300, 23700, 21300, 19000, 18000, 19000, 20400, 19200, 16200, 16100, 13400, 15400, 29500, 32800, 50500, 60400, 40400, 47400)
ABC <- c(44200, 37100, 33400, 28800, 25200, 25000, 28800, 25300, 19600, 17200, 16800, 15900, 17200, 16900, 17300, 20900, 23000, 21000, 21000, 20100, 18000, 16100, 15200, 16000, 17200, 16200, 13700, 13700, 11800, 13100, 15000, 15100, 22000, 29600, 34500, 40500)
TAC <- c(18000, 19300, 17300, 14500, 14800, 13500, 21400, 27700, 36400, 32200, 33200, 28800, 25200, 25000, 28800, 25300, 19400, 16800, 16800, 15400, 17200, 16900, 17300, 20900, 22600, 21000, 20700, 20100, 18000, 16100, 15200, 16000, 17200, 16200, 13700, 13700, 11800, 13100, 15000, 15100, 18300, 26100, 34500, 39600)
lan <- c(10400, 12600, 12000, 11800, 14100, 14500, 28900, 35200, 38400, 34800, 30200, 26400, 23900, 25400, 23600, 20700, 17400, 14600, 13900, 13600, 15600, 14100, 14700, 16400, 17500, 16600, 15600, 16000, 14600, 13100, 11900, 13000, 13900, 13600, 11500, 10900, 10200, 12300, 14200, 16600, 19000, 21300, 26900, 20400)

names(OFL) <- 1992:2023
names(ABC) <- 1988:2023
names(TAC) <- 1980:2023
names(lan) <- 1980:2023

plot(1980:2023, TAC, type="l", col="black", ylim=c(0, 75000))
lines(1988:2023, ABC, col="red")
lines(1992:2023, OFL, col="blue")
lines(1980:2023, lan, col="purple")

management_quantities <- data.frame(
    year = 1980:2023,
    OFL=c(rep(NA, length(TAC) - length(OFL)), OFL),
    ABC=c(rep(NA, length(TAC) - length(ABC)), ABC),
    TAC=c(rep(NA, length(TAC) - length(TAC)), TAC),
    landings=c(rep(NA, length(TAC) - length(lan)), lan)
) %>% pivot_longer(OFL:landings, names_to="quantity", values_to="value") %>%
    mutate(
        quantity = factor(quantity, levels=c("OFL", "ABC", "TAC", "landings"), labels=c("OFL", "ABC", "TAC", "Landings"))
    )

ggplot(management_quantities, aes(x=year, y=value, color=quantity))+
    geom_line(size=1.5)+
    scale_y_continuous(limits=c(0, 75000), labels=scales::comma)+
    labs(x="Year", y="Catch (mt)", title="Alaska Sablefish Management Quantities", color="Quantity")+
    coord_cartesian(expand=0)+
    theme_bw()


management_quantities %>%
    group_by(quantity) %>%
    mutate(
        change = abs((value-lag(value, 1))/value)
        # change = replace_na(change, 0)
    ) %>%
    filter(year > 1992) %>%
    summarise(
        trigger_10 = sum(change >= 0.10, na.rm=TRUE),
        trigger_15 = sum(change >= 0.15, na.rm=TRUE),
        trigger_25 = sum(change >= 0.25, na.rm=TRUE)
    ) %>%
    print(n=100)
