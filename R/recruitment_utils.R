#' Resample recruitment from historical timeseries
#' #'
#' Resample, with replacement, the historical recruitment timeseries.
#'
#' @param hist_recruits historical timeseries of recruitment (or deviates)
#' @param nyears total number of years to resample
#' @param seed random seed for reproducability
#'
#' @export resanple_recruits
#'
#' @example
#'
resample_recruits <- function(hist_recruits, nyears, seed){
    set.seed(seed)
    return(sample(hist_recruits, size=nyears, replace=TRUE))
}

#' Regime-like recruitment
#' #'
#' Generate future recruitment timeseries based on historcial regime-like 
#' recruitment dynamics. Historical recruitments identified as having belonged
#' to a specific regime are randomly resampled, with replacement.
#'
#' @param regime1_recruits historical recruitments from first regime
#' @param regime2_recruits historical recruitments from second regime
#' @param nyears total number of years to resample
#' @param regime_length vector inidcating how long each regime last (limited to 2 regimes)
#' @param starting_regime regime to start with (0 = first regime, 1 = second regime)
#' @param seed random seed for reproducability
#' 
#' @export resample_regime_recruits
#'
#' @example
#'
resample_regime_recruits <- function(regime1_recruits, regime2_recruits, nyears, regime_length, starting_regime, seed){
    set.seed(seed)
    rec <- rep(NA, nyears)
    curr_regime <- starting_regime
    y <- 1
    while(y < length(rec)){
        reg_len <- regime_length[curr_regime+1]
        if(curr_regime == 0){
            rs <- sample(regime1_recruits, size=reg_len, replace=TRUE)
        }else{
            rs <- sample(regime2_recruits, size=reg_len, replace=TRUE)
        }
        rec[y:(y+reg_len-1)] <- rs
        curr_regime <- !curr_regime
        y <- y+reg_len
    }
    return(rec)
}


#' Regime-like recruitment
#' 
#' Generate future recruitment timeseries based on mean and CV of distinct
#' regimes, via a lognormal distribution.
#'
#' @param mus vector of mean recruitment in each regime
#' @param cvs vector of CV of recruitment in each regime (CVs are converted internally to log SD)
#' @param nyears total number of years to resample
#' @param regime_length vector inidcating how long each regime last (limited to 2 regimes)
#' @param starting_regime regime to start with (0 = first regime, 1 = second regime)
#' @param seed random seed for reproducability
#' 
#' @export resample_regime_recruits
#'
#' @example
#'
regime_recruits <- function(mus, cvs, nyears, regime_length, starting_regime, seed){
    set.seed(seed)
    rec <- rep(NA, nyears)
    curr_regime <- starting_regime
    y <- 1
    while(y < length(rec)){
        reg_len <- regime_length[curr_regime+1]
        sds <- sqrt(log(cvs[curr_regime+1]^2 + 1))
        rec[y:(y+reg_len-1)] <- rlnorm(reg_len, meanlog=log(mus[curr_regime+1])-(sds^2)/2, sdlog=sds)
        curr_regime <- !curr_regime
        y <- y+reg_len
    }
    return(rec)
}

#' Beverton-Holt stock recruitment relationship
#' #'
#' Defines a function factory for yieding recruitment from a Bevrton-Holt
#' stock recruitment relationship, parameterized using steepeness. Initial
#' call to this function sets up a new function that takes the current SSB
#' as input, and returns the predicted recruitment from the BH SRR defined
#' by h, R0, and S0.
#'
#' @param h steepness (0.2 <= h <= 1.0)
#' @param R0 unfished recruitment
#' @param S0 unfished spawning biomass
#'
#' @export beverton_holt
#'
#' @examples
#' \dontrun{
#'  bevholt <- beverton_holt(0.7, 25, 300)
#'  bevholt(ssb=185) # 23.43891
#' }
#'
beverton_holt <- function(h, R0, S0, sigR, seed){
    # note that the set.seed() call needs to happen
    # outside of the returned function, or else there
    # will be no random variability in recruitment draws
    set.seed(seed)
    function(ssb){
        bh <- (4*R0*h*ssb)/((1-h)*R0*(S0/R0) + (5*h - 1)*ssb)
        return(
            bh + rnorm(1, mean=0, sd=sigR)
        )
    }
}

#' Recruitment function with a specified crash period
#' 
#' Generate future recruitment timeseries via historical reasampling, 
#' (see `resample_recruits`) but set a specific time period of the
#' projected recruitment to a different values. Can be used to specify
#' a "crash" in recruitment, or a "spike" in recruitment depending on
#' the `crash_value`.
#'
#' @param crash_start_years projection year to start the recruitment crash
#' @param crash_length number of years to have "crashed" recruitment
#' @param crash_value recruitment during the "crash" period
#' @param hist_recruits historical timeseries of recruitment (or deviates)
#' @param nyears total number of years to resample
#' @param seed random seed for reproducability
#'
#' @export recruits_crash
#'
#' @example
#'
recruits_crash <- function(crash_start_year, crash_length, crash_value, hist_recruits, nyears, seed){
    rec <- resample_recruits(hist_recruits, nyears, seed)
    rec[crash_start_year:(crash_start_year+crash_length-1)] <- crash_value
    return(rec)
}

retro_recruitment <- function(hist_recruits, nyears, seed, sim_func, ...){
    rec <- rep(NA, nyears)
    rec[1:length(hist_recruits)] <- hist_recruits
    rec[(length(hist_recruits)+1):nyears] <- do.call(sim_func, c(list(hist_recruits=hist_recruits), list(nyears=nyears), list(seed=seed), list(...)))#rlnorm(nyears-length(rec[!is.na(rec)]), meanlog=log(15), sdlog=0.20)
    return(rec)
}

f <- function(...){
    return(list(...))
}

f(a=1, b=2)
