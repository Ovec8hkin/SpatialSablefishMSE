simulate_weighted_comps_ISS <- function(aa_matrix, selex, region_weights, total_samples, aggregate_sex, as_integers, return_agg=TRUE){

    std_region_weights <- as.vector(region_weights/sum(region_weights))
    
    ac <- afscOM::simulate_ac(aa_matrix, selex, aggregate_sex=aggregate_sex)

    regional_samples <- round(total_samples*std_region_weights, 0)

    nregions <- afscOM::get_model_dimensions(selex)$nregions
    ac_obs_out <- array(NA, dim=c(1, 30, 2, nregions, 1))
    for(r in 1:nregions){
        tmp <- afscOM::simulate_multinomial_obs(
            ac[,,,r,drop=FALSE], 
            regional_samples[r], 
            aggregate_sex = aggregate_sex,
            as_integers=as_integers
        )
        ac_obs_out[1,,,r,1] <- tmp
    }
    
    if(return_agg){
        return(array(apply(ac_obs_out, c(1, 2, 3, 5), sum), dim=c(1,30,2,1,1)))
    }else{
        return(ac_obs_out)
    }

}


simulate_weighted_comps_SAMPLE <- function(aa_matrix, selex, region_weights, total_samples, aggregate_sex, as_integers){

    std_region_weights <- as.vector(region_weights/sum(region_weights))
    
    # aa_total <- array(apply(aa_matrix, c(1, 2, 3), sum), dim=c(1, 30, 2, 1))

    ac <- afscOM::simulate_ac(aa_matrix, selex, aggregate_sex=aggregate_sex)
    
    aa_total <- array(apply(ac*std_region_weights, c(1, 2, 3), sum), dim=c(1, 30, 2, 1))

    tmp <- afscOM::simulate_multinomial_obs(
        aa_total, 
        total_samples, 
        aggregate_sex = aggregate_sex,
        as_integers=as_integers
    )

    ac_obs_out <- array(NA, dim=c(1, 30, 2, 1, 1))
    ac_obs_out[1,,,1,1] <- tmp
    return(ac_obs_out)

}

simulate_weighted_comps_HYRBID <- function(aa_matrix, selex, region_weights, total_samples, aggregate_sex, as_integers){
    std_region_weights <- as.vector(region_weights/sum(region_weights))
    ISS_weighted_comps <- simulate_weighted_comps_ISS(aa_matrix, selex, region_weights, total_samples, aggregate_sex, as_integers, return_agg=FALSE)
    regional_samples <- round(total_samples*std_region_weights, 0)

    weighted_comp <- sweep(ISS_weighted_comps, 4, std_region_weights, FUN="*")
    std_weighted_comp <- array(NA, dim=dim(weighted_comp))
    for(r in 1:length(regional_samples)){
        tmp <- weighted_comp[,,,r,,drop=FALSE]/sum(weighted_comp[,,,r,,drop=FALSE])
        std_weighted_comp[,,,r,] <- tmp
    }

    full_weighted_comp <- sweep(std_weighted_comp, 4, regional_samples, FUN="*")
    
    total_comp <- array(apply(full_weighted_comp, c(1, 2, 3), sum), dim=c(1, 30, 2, 1))
    return(total_comp)
}
