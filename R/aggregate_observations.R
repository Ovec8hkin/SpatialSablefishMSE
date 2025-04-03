aggregate_observations <- function(survey_obs, y){

    aggregate_index <- function(index){
        return(sum(index))
    }

    aggregate_comps <- function(comps){
        return(sum(comps))
    }

    out_dims1 <- dim(survey_obs$rpns)
    out_dims1[1] <- y
    out_dims1[4] <- 1

    out_dims2 <- dim(survey_obs$acs)
    out_dims2[1] <- y
    out_dims2[4] <- 1

    aggregate_obs <- list(
        rpns = array(apply(survey_obs$rpns, c(1, 5), \(x) aggregate_index(x))[1:y,], dim=out_dims1),
        rpws = array(apply(survey_obs$rpws, c(1, 5), \(x) aggregate_index(x))[1:y,], dim=out_dims1),
        acs = array(apply(survey_obs$acs, c(1, 2, 3, 5), \(x) aggregate_comps(x))[1:y,,,], dim=out_dims2)
    )

    return(aggregate_obs)

}
