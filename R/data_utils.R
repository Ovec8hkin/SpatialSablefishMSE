#' Change a `ggdist` Tibble to Long Format 
#' #'
#' Description
#'
#' @param data a tibble output by a `ggdist` summary function
#' (e.g. `ggdist:median_qi(...)`)
#' @param n the number of grouping variables used (default: 1)
#'
#' @return a long format tibble with columns median, lower, and upper
#' @export reformat_ggdist_long
#'
#' @example
#'
reformat_ggdist_long <- function(data, n=1){

    data_long <- data %>% pivot_longer(-c(1:n, ncol(.), ncol(.)-1, ncol(.)-2))
    medians <- data_long %>% filter(!(grepl("lower",name) | grepl("upper",name)))
    lowers <- data_long %>% filter(grepl("lower", name)) %>% mutate(name=unlist(lapply(str_split(name, fixed(".")), `[[`, 1)))
    uppers <- data_long %>% filter(grepl("upper", name)) %>% mutate(name=unlist(lapply(str_split(name, fixed(".")), `[[`, 1)))

    return(
        medians %>% 
            left_join(lowers, by=c(colnames(medians)[1:(5+n-1)])) %>% 
            left_join(uppers, by=c(colnames(medians)[1:(5+n-1)])) %>%
            rename("median"="value.x", "lower"="value.y", "upper"="value")
    )
    
}

#' Create Tibble of MSE Output Data from Multiple Models
#'
#' @param model_runs a list of MSE model run objects (created via `run_mse(...)`)
#' @param var the output variable from the model run objects
#' @param extra_columns a list of extra column names and values to add to the tibble
#'
#' @export bind_mse_outputs
#'
#' @example
#'
bind_mse_outputs <- function(model_runs, var, extra_columns){

    bind_rows(
        lapply(
            seq_along(model_runs), 
            function(x){
                y <- melt(model_runs[[x]][var])
                for(i in 1:length(extra_columns)){
                    y <- y %>% mutate(!!names(extra_columns)[[i]] := extra_columns[[i]][x])
                }
                return(y)
            }
        )
    )

}