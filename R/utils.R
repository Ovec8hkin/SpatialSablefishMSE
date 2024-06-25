#' Calculate Joint Selectivity and Retention Across Multiple Fleets
#' 
#' Computes the average selectivity-at-age and retention-at-age acting
#' on a population when multiple fleets are present. Selectivity and
#' retention are weighted based on user supplied weights.
#'
#' @param sel selectiviity-at-age ([1, nages, nesexes, nregions, nfleets])
#' @param ret retention-at-age ([1, nages, nsexes, nregions, nfleets])
#' @param prop_fs fleet weights
#'
#' @export calculate_joint_selret
#'
#' @example
#'
calculate_joint_selret <- function(sel, ret, prop_fs=c(0.50, 0.50)){
    joint_self <- apply(sel[,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum)/max(apply(sel[,,1,,,drop=FALSE]*prop_fs, c(1, 2), sum))
    joint_selm <- apply(sel[,,2,,,drop=FALSE]*prop_fs, c(1, 2), sum)/max(apply(sel[,,2,,,drop=FALSE]*prop_fs, c(1, 2), sum))
    joint_retf <- apply(ret[,,1,,,drop=FALSE], c(1, 2), sum)/max(apply(ret[,,1,,,drop=FALSE], c(1, 2), sum))
    joint_retm <- apply(ret[,,2,,,drop=FALSE], c(1, 2), sum)/max(apply(ret[,,2,,,drop=FALSE], c(1, 2), sum))
    
    joint_sel <- array(NA, dim=dim(sel)[1:4])
    joint_sel[,,1,] <- joint_self
    joint_sel[,,2,] <- joint_selm

    joint_ret <- array(NA, dim=dim(ret)[1:4])
    joint_ret[,,1,] <- joint_retf
    joint_ret[,,2,] <- joint_retm

    return(list(sel=joint_sel, ret=joint_ret))
}


#' Create dplyr Filter for Years Around Large Recruitment Events
#' 
#' Given a list of MSE model runs, identifies years in which recruitment was
#' larger than a specified threshold value for each OM simulation. That information
#' is used to generate an complex boolean expression object for use with 
#' dplyr::filter(...) when computing performance metrics.
#'
#' @param om_list a list of OM objects (as would be given to `run_mse_multiple`)
#' @param seed_list a vector or list of simulation seeds (as would be given to `run_mse_multiple`)
#' @param model_runs list of completed MSE model runs (the output of `run_mse_multiple`)
#' @param large_event_thresh threshold value for declaring a recruitment event "large"
#'
#' @export 
#'
#' @example
#'
get_big_recruitment_filter <- function(om_list, seed_list, model_runs, large_event_thresh){
    rec_years_interest <- list()

    for(o in 1:length(om_list)){
        oname <- names(om_list)[o]
        mr1 <- model_runs[[o]]
        for(s in 1:length(seed_list)){
            sim <- seed_list[s]
            mr1_recs <- apply(mr1$naa[65:110,1,,,s], 1, sum)
            # large_event_thresh <- 60
            large_event_years <- as.numeric(names(which(mr1_recs >= large_event_thresh)))
            years <- unique(as.vector(sapply(large_event_years, \(x) return(seq(x-2, x+5, 1)))))
            exp <- 
                paste0(
                    "om == '", oname ,
                    "' & sim == ", sim, 
                    " & time %in% c(", paste(years, collapse=", "), ")"
                )
            rec_years_interest <- c(exp, rec_years_interest)
        }
    }

    eval(rec_years_interest[1])

    giant_filter <- rlang::parse_expr(paste(paste0("(", rec_years_interest, ")"), collapse = " | "))
    return(giant_filter)
}

#' Average Annual Variation (AAV)
#' 
#' Calculates average annual variation ($\frac{\sum{\frac{|x_y - x_{y-1}|}{mean(x)}}}{N-1}$)
#'
#' @param data ordered vector of observations of a quantity
#'
#' @export aav
#'
#' @example
#'
aav <- function(data){
    total <- mean(data)
    diffs <- abs(diff(data))
    aav <- sum(diffs/total)/(length(data)-1)
    return(ifelse(is.nan(aav), 0, aav)) # If all data is 0, return 0 rather than NA
}