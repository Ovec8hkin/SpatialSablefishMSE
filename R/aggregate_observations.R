aggregate_indices <- function(survey_obs, y){

    out_dims1 <- dim(survey_obs$rpns)
    out_dims1[1] <- y
    out_dims1[4] <- 1

    aggregate_obs <- list(
        rpns = array(apply(survey_obs$rpns, c(1, 5), \(x) sum(x))[1:y,], dim=out_dims1),
        rpws = array(apply(survey_obs$rpws, c(1, 5), \(x) sum(x))[1:y,], dim=out_dims1)
    )

    return(aggregate_obs)

}

