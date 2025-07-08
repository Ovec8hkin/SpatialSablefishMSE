#' Regional ABC Apportionment by RPW Moving Average 
#' 
#' The operational method for regional ABC apportionment for sablefish
#' in Alaska is a simple 5-yr moving average across the Alaska LL survey
#' RPW index.
#' 
#' @param survey_obs list of survey observations with an RPW list element
#' @param y the terminal year of observations to calculate the moving average over
#' @param window_size size of the moving average window to use
#' 
#' @returns proportion of survey RPW in each spatial region in year y
#' 
#' @export rpw_moving_average
#' 
rpw_moving_average <- function(survey_obs, y, window_size=5){
    rpw <- survey_obs$rpws[1:y,,,,1]
    roll <- apply(rpw, 2, \(x) zoo::rollmean(x, k=window_size))
    app <- t(apply(roll, 1, \(x) x/sum(x)))
    return(app[nrow(app),])
}

#' Fixed Regional ABC Apportionment
#' 
#' Set a temporally fixed regional ABC apportionment.
#' 
#' @param app fixed regional apportionment (dim [1, nregions])
#' @param ... extra posisble parameters for generalized compatability
#' with other apportionment functions
#' 
#' @return vector of regional apportionment
#' 
#' @export fixed_apportionment
fixed_apportionment <- function(app, ...){
    return(app)
}
