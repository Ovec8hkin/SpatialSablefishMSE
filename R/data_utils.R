#' Change a `ggdist` Tibble to Long Format 
#' 
#' Pivot a summary tibble from `ggdist` with multiple summary variables
#' to long format. Summary variable names will be be pivoted to a 
#' 'name' column, and corresponding 'median', 'lower', 'upper' columns
#' will contain computed quantile values.
#' 
#' Intended for use when plotting distributions of multiple variables
#' simultaneously.
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

    if(n==0){
        return(data)
    }

    data_long <- data %>% pivot_longer(-c(1:n, ncol(.), ncol(.)-1, ncol(.)-2))
    medians <- data_long %>% filter(!(grepl("lower",name) | grepl("upper",name)))
    lowers <- data_long %>% filter(grepl("lower", name)) %>% mutate(name=unlist(lapply(str_split(name, fixed(".")), `[[`, 1)))
    uppers <- data_long %>% filter(grepl("upper", name)) %>% mutate(name=unlist(lapply(str_split(name, fixed(".")), `[[`, 1)))

    return(
        medians %>% 
            left_join(lowers, by=c(colnames(medians)[1:(4+n)])) %>% 
            left_join(uppers, by=c(colnames(medians)[1:(4+n)])) %>%
            rename("median"="value.x", "lower"="value.y", "upper"="value")
    )
    
}

#' Create Tibble of MSE Output Data from Multiple Models
#'
#' @param model_runs a list of MSE model run objects (created via `run_mse(...)`)
#' @param var the output variable from the model run objects
#' @param extra_columns a data.frame of extra column names and values to add to the tibble
#'
#' @export bind_mse_outputs
#'
#' @example
#' 
bind_mse_outputs <- function(model_runs, var, extra_columns){

    # if(length(extra_columns) == 1){
    #     col_name <- names(extra_columns)
    #     model_grid <- as.data.frame(tibble(!!col_name := extra_columns[[1]]))
    # }else{
    #     model_grid <- expand.grid(extra_columns)
    # }
    model_grid <- extra_columns
    t <- bind_rows(
        lapply(
            seq_along(model_runs), 
            function(x){
                y <- melt(model_runs[[x]][var])
                for(i in 1:ncol(model_grid)){
                    col_name <- names(model_grid)[i]
                    y <- y %>% mutate(!!col_name := model_grid[x,i])

                }
                return(y)
            }
        )
    )

}


relativize_performance <- function(data, rel_column, value_column, rel_value, grouping){
    if(is.null(rel_value)){
        return(data)
    }
    
    return(
        data %>%
            group_by(across(all_of(grouping))) %>%
            pivot_wider(names_from=rel_column, values_from = value_column) %>%
            mutate(across(everything(), ~ . / eval(rlang::parse_expr(rel_value)))) %>%
            pivot_longer((length(grouping)+1):(ncol(.)), names_to=rel_column, values_to=value_column)
    )
}

#' Filter processed MSE data to between two time points,
#' 
#' Filter dataframe of processed MSE data by 'time'
#' column to between times set by `time_horizon`.
#' 
#' @param data dataframe of processed MSE data
#' @param time_horizon vector of times to filter between.
#' If first element is NA, lower bound will be 1. If
#' second element is NA, upper bound will be maximum value
#' in `time` column of `data`.
#' 
#' @export filter_times
#' 
filter_times <- function(data, time_horizon){
    times <- seq(
        from = ifelse(!is.na(time_horizon[1]), time_horizon[1], 1),
        to   = ifelse(!is.na(time_horizon[2]), time_horizon[2], max(data$time)),
        by=1
    )

    return(data %>% filter(time %in% times))
}

