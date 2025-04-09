simulate_weighted_comps <- function(aa_matrix, selex, region_weights, total_samples, aggregate_sex, as_integers){

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
    
    return(ac_obs_out)

}
