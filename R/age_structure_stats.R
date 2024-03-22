#' Compute Shannon Diversity Index (H')
#' 
#' Calculate the Shannon diversity index (H')
#' (Shannon and Weaver, 1949) give a numbers
#' at age vector.
#' 
#' H' = sum(naa*log(naa))
#'
#' @param naa numbers-at-age vector
#'
#' @export shannon_diversity
#'
#' @example
#'
shannon_diversity <- function(naa){
    if(sum(naa) > 1) naa <- naa/sum(naa)
    return(-sum(naa*log(naa)))
}

#' Compute average age of a population
#' #'
#' Calculate the average age of a population
#' given its age structure.
#'
#' @param naa numbers-at-age vector
#' @param ages vector of ages
#'
#' @export average_age
#'
#' @example
#'
average_age <- function(naa, ages){
    return(weighted.mean(ages, naa))
}

#' Compute proportion of a population that is mature
#' 
#' Calculate the proportion of a population that
#' has reached maturity given its age structure
#' and maturity schedule.
#'
#' @param naa numbers-at-age vector
#' @param mat maturity-at-age vector
#'
#' @export prop_mature
#'
#' @example
#'
prop_mature <- function(naa, mat){
    #if(sum(naa) > 1) naa <- naa/sum(naa)
    return(sum(naa*mat)/sum(naa))
}

#' Compute proportion of population that is fully mature
#' #'
#' Calculate the proportion of a population
#' in age classes that are "fully mature" (e.g.
#' mat > 0.995).
#'
#' @param naa numbers-at-age vector
#' @param mat maturity-at-age vector
#'
#' @export prop_fully_mature
#'
#' @example
#'
prop_fully_mature <- function(naa, mat){
    if(sum(naa) > 1) naa <- naa/sum(naa)
    full_maturity <- as.numeric(mat > 0.995)
    return(sum(naa * full_maturity))
}

#' Compute Age-Based Indicator
#' #'
#' Compute ABI (Griffiths et al. 2023) given a
#' population age structure and an age structure
#' to use as a reference.
#'
#' @param naa numbers-at-age vector
#' @param ref a reference naa vector
#' @param threshold proportion of population to consider "old"
#' @param start_age age at which to start calculating ABI
#'
#' @export abi
#'
#' @example
#'
abi <- function(naa, ref, threshold=0.90, start_age=3){
    ref_naa <- ref[(start_age-1):length(ref)]
    ref_naa_prop <- ref_naa/sum(ref_naa)
    A_ref <- min(which(cumsum(ref_naa_prop) > threshold))-1
    P_ref <- sum(ref_naa_prop[(A_ref+1):length(ref_naa_prop)])

    naa <- naa[(start_age-1):length(naa)]
    naa_ref <- naa[(A_ref+1):length(naa)]
    return((sum(naa_ref)/sum(naa))/P_ref)
}
