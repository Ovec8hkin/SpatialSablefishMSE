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

#' Filter processed MSE data to a subset of HCRs and OMs
#' 
#' Filter dataframe of processed MSE data to only include
#' a subset of HCRs and OMs.
#' 
#' @param data dataframe of processed MSE data
#' @param hcrs vector of HCR names to retain in the data
#' @param oms vector of OM names to retain in the data
#' 
#' @export filter_hcr_om 

filter_hcr_om <- function(data, hcrs, oms){
    return(
        data %>% filter(hcr %in% hcrs, om %in% oms)
    )
}

#' Round dataframe column to Zero
#' 
#' Rounds data in a dataframe column to 0 when below a 
#' certain threshold. Also rounds non finite values to 
#' 0 automatically.
#' 
#' @param data dataframe of processed MSE data
#' @param col_name name of data column to round
#' @param zero_threshold threshold to round down to zero
#' 
#' @export round_to_zero
#' 
round_to_zero <- function(data, col_name, zero_threshold=1e-2){
    return(
        data %>%
            mutate(!!col_name := 
                ifelse(eval(rlang::parse_expr(col_name)) < zero_threshold, 0, eval(rlang::parse_expr(col_name)))
            ) %>%
            mutate(!!col_name := 
                ifelse(!is.finite(eval(rlang::parse_expr(col_name))), 0, eval(rlang::parse_expr(col_name)))
            )
    )
}

#' Scale and Rank dataframe columns within groups
#' 
#' Computes 0-1 scaled version of dataframe column, within
#' groups, and ranks (1-N) entries within group based on
#' scaled values (where 1 is highest in-group scaled value, and
#' N is the lowest in-group scaled value). 
#' 
#' Intended for use on performance metric data only. Note that
#' "Catch AAV, "Proportion of Years with Low SSB", and "Recovery
#' Time" performance metrics use inverse ranking system (smaller
#' scaled values are better), as the intention is to minimize
#' these metrics.
#' 
#' @param data dataframe of processed MSE data
#' @param col_name name of data column to scale and rank
#' 
#' @export scale_and_rank
#' 
scale_and_rank <- function(data, col_name){
    return(
        data %>%
            mutate(
                scaled = eval(rlang::parse_expr(col_name))/inf_max(eval(rlang::parse_expr(col_name)))
            ) %>%
            arrange(desc(scaled), .by_group=TRUE) %>%
            mutate(
                rank = ifelse(
                            name %in% c(
                                "Catch AAV", 
                                "Proportion of Years with Low SSB", 
                                "Recovery Time"
                            ), 
                            factor(desc(row_number())), 
                            factor(row_number())
                        ),
                rank = factor(rank)
            )
    )
}

