npfmc_tier3_F <- function(ssb, b40, f40, alpha){
    if(ssb/b40 >= 1)  F <- f40 # if stock status >= 1
    if(ssb/b40 > alpha && ssb/b40 < 1) F <- f40 * ((ssb/b40)-alpha)/(1-alpha) # if stock status > alpha & stock status < 1
    if(ssb/b40 <= alpha) F <- 0 # if stock stats <= alpha
    
    return(F)
}

# npfmc_tier3a_F(ssb=10, fmax=F40)
