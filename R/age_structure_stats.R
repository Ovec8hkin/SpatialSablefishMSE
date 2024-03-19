shannon_diversity <- function(naa){
    if(sum(naa) > 1) naa <- naa/sum(naa)
    return(-sum(naa*log(naa)))
}

average_age <- function(naa, ages){
    return(weighted.mean(ages, naa))
}

prop_mature <- function(naa, mat){
    #if(sum(naa) > 1) naa <- naa/sum(naa)
    return(sum(naa*mat)/sum(naa))
}

prop_fully_mature <- function(naa, mat){
    if(sum(naa) > 1) naa <- naa/sum(naa)
    full_maturity <- as.numeric(mat > 0.995)
    return(sum(naa * full_maturity))
}

abi <- function(naa, ref_naa, threshold=0.90){
    ref_naa <- ref_naa[2:length(ref_naa)]
    ref_naa_prop <- ref_naa/sum(ref_naa)
    A_ref <- min(which(cumsum(ref_naa_prop) > threshold))-1
    P_ref <- sum(ref_naa_prop[(A_ref+1):length(ref_naa_prop)])

    naa <- naa[2:length(naa)]
    naa_ref <- naa[(A_ref+1):length(naa)]
    return((sum(naa_ref)/sum(naa))/P_ref)
}
