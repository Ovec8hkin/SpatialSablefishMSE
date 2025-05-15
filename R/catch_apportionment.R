rpw_moving_average <- function(survey_obs, y, window_size=5){
    rpw <- survey_obs$rpws[1:y,,,,1]
    roll <- apply(rpw, 2, \(x) zoo::rollmean(x, k=window_size))
    app <- t(apply(roll, 1, \(x) x/sum(x)))
    return(app[nrow(app),])
}
