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
#' #'
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
