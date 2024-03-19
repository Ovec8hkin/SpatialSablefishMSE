npfmc_tier3a_F <- function(ssb, fmax){
    lrp <- 0.05
    urp <- 1.0
    F <- 0.0
    if(ssb >= urp){
        F <- fmax
    }else if(lrp < ssb && ssb < urp){
        F <- fmax * (ssb-lrp)/(urp-lrp)
    }
    return(F)
}

constant_F <- function(F){
    return(F)
}
# npfmc_tier3a_F(ssb=10, fmax=F40)
