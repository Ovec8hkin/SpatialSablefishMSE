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
    joint_self <- apply(sweep(sel[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(sel[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    joint_selm <- apply(sweep(sel[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(sel[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    joint_retf <- apply(sweep(ret[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(ret[,,1,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    joint_retm <- apply(sweep(ret[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(ret[,,2,,,drop=FALSE], 5, prop_fs, FUN="*"), c(1, 2), sum))
    
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
get_big_recruitment_filter <- function(om_list, seed_list, model_runs, om_names, large_event_thresh, lags){
    rec_years_interest <- list()

    for(o in 1:length(om_list)){
        oname <- om_names[o]
        mr1 <- model_runs[[o]]
        for(s in 1:length(seed_list)){
            sim <- seed_list[s]
            mr1_recs <- apply(mr1$naa[65:110,1,,,s], 1, sum)
            # large_event_thresh <- 60
            large_event_years <- as.numeric(names(which(mr1_recs >= large_event_thresh)))
            years <- unique(as.vector(sapply(large_event_years, \(x) return(seq(x+lags[1], x+lags[2], 1)))))
            exp <- 
                paste0(
                    "om == '", oname ,
                    "' & sim == ", sim, 
                    " & time %in% c(", paste(years, collapse=", "), ")"
                )
            rec_years_interest <- c(exp, rec_years_interest)
        }
    }

    # eval(rec_years_interest[1])

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

#' Extend a 3D array along its final dimension
#' 
#' Copies and append values from the final dimension 
#' of a 3D array to the end of the final dimension.
#' Originally written by Craig Marsh.
#' 
#' @param array_3d array with 3 dimensions
#' @param n number of times to concatenate the last element in 
#' the 3rd dimension
#' 
#' @export extend_3darray_last_dim
#' 
extend_3darray_last_dim <- function(array_3d, n){
    if(n <= 0){
        return(array_3d[,,1:(dim(array_3d)[3]+n)])
    }
    new_3d_array = abind(array_3d, replicate(array_3d[, , dim(array_3d)[3]], 
        n = n), along = 3)
    return(new_3d_array)
}

#' Extend a vector by copying its final value
#' 
#' Copies and appends the final value of a vector
#' to the end of said vector 'n' times.
#' 
#' @param vector a vector to add n elements to
#' @param n number of times to concatenate final element
#' 
#' @export extend_vec_last_val
#' 
extend_vec_last_val <- function(vector, n){
    if(n <= 0){
        return(vector[1:(length(vector)+n)])
    }
    new_vec = c(vector, rep(vector[length(vector)], n))
    return(new_vec)
}

#' Generate frequency vector
#' 
#' Create a vector of 0s and 1s that that alternate with
#' a given frequency. This is intended to be used to
#' indicate when certain actions should occur within a 
#' loop.
#' 
#' @param frequency freqency with which a 1 occurs
#' @param len length of output vector
#' 
#' @export generate_annual_frequency
#' 
generate_annual_frequency <- function(frequency, len){
    do <- rep(1, len+1)
    if(length(frequency) > 1){
        do <- frequency
    }else if(length(frequency) == 1){
        survey_years <- rep(0, len+1)
        survey_years[seq(1, length(survey_years), frequency)] <- 1
        survey_years[length(survey_years)+1] <- 1
        do <- survey_years
    }
    return(do)
}

#' Load Saved MSE Model Runs from Disk
#' 
#' Read all saved RDS files present in data/active and coerce into proper
#' model_runs list object. Also setup correctly specified extra_columns
#' object for use with bind_mse_outputs. 
#'
#' @param om_order vector of correct order of OMs (used to set OM factor level)
#' @param hcr_order vector of correct order of HCRs (used to set HCR factor level)
#'
#' @export get_saved_model_runs
#'
#' @example
#'
get_saved_model_runs <- function(om_order=NULL, hcr_order=NULL){
    fs <- list.files(file.path(here::here(), "data", "active"), full.names = TRUE)
    model_runs <- unlist(lapply(fs, function(x){
        m <- readRDS(x)
        mse <- m$mse_objects
        mse[(length(mse)-3):length(mse)]
    }), recursive=FALSE)

    om_names <- sapply(fs, function(x){
        m <- readRDS(x)
        lapply(m$om_list, function(om){
            om$name
        })
    })
    om_names_formatted <- unlist(c(om_names))

    hcr_names <- sapply(fs, function(x){
        m <- readRDS(x)
        m$hcr$name
    })
    hcr_names_repped <- rep(hcr_names, each=length(unique(om_names_formatted)))

    extra_columns2 <- data.frame(om=om_names_formatted, hcr=hcr_names_repped)

    if(!is.null(om_order))
        extra_columns2$om <- factor(extra_columns2$om, levels=om_order)
    
    if(!is.null(hcr_order))
        extra_columns2$hcr <- factor(extra_columns2$hcr, levels=hcr_order)

    return(listN(model_runs, extra_columns2))

}

#' Get maximum value without considering infinite values
#' 
#' Wrapper around max that ignores infinite values
#' 
#' @param d vector of values to find maximum of
#' 
#' @export inf_max
#' 
inf_max <- function(d){
    return(max(d[!is.infinite(d)]))
}
