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

    sweep_dim <- 5
    if(length(dim(prop_fs))==2){
        sweep_dim <- c(1, 5)
    }else if(length(dim(prop_fs))==3){
        sweep_dim <- c(1, 4, 5)
    }else{
        sweep_dim <- 5
    }
        
    if(!all(prop_fs == 0)){
        joint_self <- apply(sweep(sel[,,1,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(sel[,,1,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum))
        joint_selm <- apply(sweep(sel[,,2,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(sel[,,2,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum))
        joint_retf <- apply(sweep(ret[,,1,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(ret[,,1,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum))
        joint_retm <- apply(sweep(ret[,,2,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum)/max(apply(sweep(ret[,,2,,,drop=FALSE], sweep_dim, prop_fs, FUN="*"), c(1, 2), sum))
    }else{
        joint_self <- array(0, dim=dim(sel)[1:2])
        joint_selm <- joint_self
        joint_retf <- array(0, dim=dim(ret)[1:2])
        joint_retm <- joint_retf
        
    }
    
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
get_saved_model_runs_old <- function(om_order=NULL, hcr_order=NULL){
    fs <- list.files(file.path(here::here(), "data", "mse_runs_good"), full.names = TRUE)
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
    if(!is.null(hcr_order))
        fs <- unlist(sapply(hcr_order, \(x) fs[grepl(sub("/", "", sub(" ", "_", tolower(x))), fs)]))
    
    model_runs <- lapply(seq_along(fs), function(i){
        x <- fs[i]
        m <- readRDS(x)
        mse <- m$mse_objects
        mse[[length(mse)]]
    })

    om_names <- sapply(fs, function(x){
        m <- readRDS(x)
        m$om$name
    })
    names(om_names) <- NULL

    hcr_names <- sapply(fs, function(x){
        m <- readRDS(x)
        m$hcr$name
    })
    names(hcr_names) <- NULL

    extra_columns2 <- data.frame(om=om_names, hcr=hcr_names)

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

#' Get minimum value without considering infinite values
#' 
#' Wrapper around max that ignores infinite values
#' 
#' @param d vector of values to find maximum of
#' 
#' @export inf_min
#' 
inf_min <- function(d){
    return(min(d[!is.infinite(d)]))
}

#' Realized Sample Size
#' 
#' Calculate realized sample size from Hulson and Williams, 2024
#' 
#' @param obs age composition observations
#' @param exp expected age composition
#' 
#' @return realized sample size sum(obs*(1-obs))/sum((exp-obs)^2)
#' 
rss <- function(obs, exp){
    num <- sum(obs*(1-obs))
    denom <- sum((exp-obs)^2)
    return(num/denom)
}

#' Realized Sample Size
#' 
#' Wrapper function for calculating realized sample size from Hulson
#' and williams, 2024.
#' 
#' @param ac_obs age composition observations
#' @param naa population age structure
#' @param selex selectivity-at-age
#' 
#' @return vector of realized sample sizes for each year
#' 
calculate_realized_samplesize <- function(ac_obs, naa, selex){
    nyears <- dim(selex)[1]
    nfleets <- dim(selex)[5]
    
    selected2 <- naa*selex
    total_selected2 <- array(apply(selected2, c(1, 2, 3), sum), dim=c(nyears+1,30,2,1))
    selected_prop2 <- aperm(apply(total_selected2, c(1), \(x) x/sum(x)), c(2, 1))

    ac2 <- array(ac_obs[,,,1,1], c(nyears, 30*2, 1))
    ac2 <- aperm(apply(ac2, c(1), \(x) x/sum(x)), c(2, 1))


    r_sampsize <- sapply(1:nyears, function(x){
        rss(obs=selected_prop2[x,], exp=ac2[x,])
    })

    r_ISS <- round(r_sampsize, 0)
    return(r_ISS)
}

#' Create Time-Varying Movement Matrix for Climate Forced OMs
#' 
#' Generate a time-varying movement matrix that linearly
#' scales between a starting and ending movement matrix
#' over a given period time. Starting matrix is from 
#' Cheng et al. 2025, while ending matrix was designed to
#' increase westward movement of adult/old fish.
#' 
#' @param time_trend duration over which to move between matrices
#' @param nyears total number of years in output movement matrix
create_tv_movement <- function(time_trend=30, nyears=200){

    start_movement_young <- matrix(
        c(
            0.792, 0.033, 0.047, 0.089, 0.039, 
            0.068, 0.820, 0.030, 0.058, 0.024, 
            0.047, 0.066, 0.665, 0.160, 0.062, 
            0.025, 0.025, 0.039, 0.829, 0.082, 
            0.017, 0.042, 0.057, 0.214, 0.670
        ),
        nrow=5, byrow=TRUE
    )

    end_movement_young <- start_movement_young

    start_movement_adu <- matrix(
        c(
            0.548, 0.059, 0.100, 0.167, 0.122, 
            0.037, 0.705, 0.049, 0.085, 0.122, 
            0.039, 0.045, 0.514, 0.264, 0.167, 
            0.012, 0.014, 0.026, 0.807, 0.139, 
            0.014, 0.023, 0.022, 0.167, 0.772
        ),
        nrow=5, byrow=TRUE
    )

    end_movement_adu <- matrix(
        c(
            0.8500, 0.0198, 0.0335, 0.0559, 0.0408, 
            0.0253, 0.8000, 0.0334, 0.0580, 0.0833, 
            0.0427, 0.1062, 0.6000, 0.1538, 0.0973, 
            0.0837, 0.0965, 0.1308, 0.5500, 0.1390, 
            0.0817, 0.0985, 0.1133, 0.2095, 0.5000
        ),
        nrow=5, byrow=TRUE
    )

    start_movement_old <- matrix(
        c(
            0.559, 0.079, 0.120, 0.124, 0.117, 
            0.041, 0.771, 0.036, 0.063, 0.087, 
            0.043, 0.049, 0.517, 0.220, 0.169, 
            0.010, 0.018, 0.023, 0.796, 0.145, 
            0.018, 0.030, 0.030, 0.143, 0.786
        ),
        nrow=5, byrow=TRUE
    )

    end_movement_old <- matrix(
        c(
            0.8500, 0.0469, 0.0409, 0.0323, 0.0299, 
            0.0361, 0.8000, 0.0317, 0.0555, 0.0767, 
            0.0468, 0.1106, 0.6000, 0.1372, 0.1054, 
            0.0798, 0.1139, 0.1177, 0.5500, 0.1386,
            0.0785, 0.0936, 0.1149, 0.2130, 0.5000
        ),
        nrow=5, byrow=TRUE
    )

    time_trend <- 30

    move_matrix <- array(NA, dim=c(5, 5, nyears, 30, 2))
    move_matrix[,,1,1:6,] <- start_movement_young
    move_matrix[,,1,7:15,] <- start_movement_adu
    move_matrix[,,1,16:30,] <- start_movement_old

    move_matrix[,,(1+time_trend):nyears,1:6,] <- end_movement_young
    move_matrix[,,(1+time_trend):nyears,6:15,] <- end_movement_adu
    move_matrix[,,(1+time_trend):nyears,16:30,] <- end_movement_old

    total_start_move <- move_matrix[,,1,,]
    total_end_move <- move_matrix[,,1+time_trend,,]

    for(i in 1:time_trend){
        change_matrix <- (total_end_move - total_start_move)/time_trend
        
        move_matrix[,,1+i,,] <- move_matrix[,,1+(i-1),,] + change_matrix
    }

    return(move_matrix)
}

#' Check MLE Convergence via Post-Optimization Sanity Checks
#' 
#' @param sd_rep standard deviation report from TMB
#' @param rep report object from TMB
#' @param gradient_tol tolerance for maximum absolute gradient
#' @param se_tol tolerance for maximum standard error
#' @param corr_tol tolerance for maximum correlation between parameters
#' @return TRUE if all checks are passed, FALSE otherwise
#' @export post_optim_sanity_checks
#'
post_optim_sanity_checks <- function(sd_rep, rep, gradient_tol = 1e-3, se_tol = 100, corr_tol = 0.99) {
    passed_post_sanity_checks <- TRUE
    
    # check likelihoods are all finite and not NA
    if(!all(is.finite(rep$jnLL))) {
        message("Found Inf in joint log-likelihood, model is not converged!")
        passed_post_sanity_checks <- F
    }
    # check maximum absolute gradients
    max_abs_grad_ndx <- which.max(abs(sd_rep$gradient.fixed))
    max_abs_grad <- abs(sd_rep$gradient.fixed)[max_abs_grad_ndx]
    if(gradient_tol < max_abs_grad) {
        message("Parameter: ", names(sd_rep$par.fixed)[max_abs_grad_ndx], " had absolute gradient = ", max_abs_grad,
                " which was greater than tolerance ", gradient_tol,". This indicates potential non-convergence according to the tolerance.\n")
        passed_post_sanity_checks <- F
    }
    # check hessian
    if(!sd_rep$pdHess) {
        message("Hessian is not positive definite, model is not converged!")
        passed_post_sanity_checks <- F
    }
    # check if standard errors are finite
    if(!all(is.finite(sqrt(diag(sd_rep$cov.fixed))))) {
        message("Found non finite elements in standard errors of parameters, model is not converged!")
        passed_post_sanity_checks <- F
    }
    # check if standard errors are big
    if(max(sqrt(diag(sd_rep$cov.fixed))) > se_tol) {
        message("Parameter: ", names(diag(sd_rep$cov.fixed))[which.max(sqrt(diag(sd_rep$cov.fixed)))], " has a standard error = ",
                max(sqrt(diag(sd_rep$cov.fixed))), " which was greated than tolerance ", se_tol, ". This indicates potential non-convergence according to the tolerance. \n")
        passed_post_sanity_checks <- F
    }
    # check if correlations are big
    corr_mat <- cov2cor(sd_rep$cov.fixed)
    diag(corr_mat) <- "Same" # set diagonal to "Same" to remove from max calculations
    # reshape to dataframe
    corr_df <- reshape2::melt(corr_mat) %>%
        dplyr::filter(value != 'Same') %>%
        dplyr::mutate(value = as.numeric(value))
    if(max(abs(corr_df$value)) > corr_tol) {
        message("Parameter pairs: ", corr_df$Var1[which.max(abs(corr_df$value))], " and ", corr_df$Var2[which.max(abs(corr_df$value))], " have a correlation of ", max(abs(corr_df$value)), ". This indicates potential non-convergence according to the tolerance.")
        passed_post_sanity_checks <- F
    }
    # cat("\n\n");
    if(passed_post_sanity_checks) {
        message("Successfully passed post-optim-sanity checks\n")
    }

    return(passed_post_sanity_checks)

}

#' Plot Diagnostics for a Given Simulation from Saved MSE Runs
#' 
#' @param hcr_filter vector of HCR names to filter on (e.g., c("HCR1", "HCR2"))
#' @param om_filter vector of OM names to filter on (e.g., c("OM1", "OM2"))
#' @param simulation_seed seed of simulation to plot
#' @param simulation_year year of simulation to plot
#' @export check_diagnostics
#'
# check_diagnostics <- function(hcr_filter, om_filter, simulation_seed, simulation_year){
#     x <- find_model_run(hcr_filter, om_filter, simulation_seed)
#     m <- x$m
#     model_run <- x$model_run
    
#     n_proj_years <- m$mse_options_list$mse_options$n_proj_years
#     n_spinup_years <- m$mse_options_list$mse_options$n_spinup_years
#     seeds <- m$seeds
#     nsims <- length(seeds)
#     simulation_number <- which(seeds == simulation_seed)
#     p <- plot_diags(model_run, n_spinup_years, n_proj_years, simulation_year=simulation_year, simulation_number = simulation_number)
#     return(p)
# }

#' Get Unconverged Model Diagnostics
#' 
#' Get parameter and problem information for unconverged models.
#' @param convergence_table data frame from unconverged models, from `check_em_convergence_diagnostics`
#' 
#' @return data frame with parameter and problem columns added
#' @export get_unconverged_model_diags
#'
get_unconverged_model_diags <- function(convergence_table, om_list, seed_list){
    convergence_table <- convergence_table %>% mutate(parameter=NA, problem=NA)
    for(i in 1:nrow(convergence_table)){
        print(paste0(i,"/", nrow(convergence_table)))
        r <- convergence_table[i, ]
        hcr_filter <- r$hcr
        om_filter <- r$om
        simulation_seed <- r$sim
        simulation_year <- r$sim_year

        x <- find_model_run(hcr_filter, om_filter, simulation_seed, om_list, seed_list)
        model_run <- x$model_run

        seeds <- x$m$seeds
        nsims <- length(seeds)
        simulation_number <- which(simulation_seed==seeds)
        n_proj_years <- x$m$mse_options_list$mse_options$n_proj_years

        em_model_obj <- model_run$model_outs[[(simulation_year)+(simulation_number-1)*(n_proj_years+1)]]

        report <- em_model_obj$rep

        # convergence_message = capture.output(type=c("message"),{
        #     sdrep = TMB::sdreport(em_model_obj)
        #     SPoRC::post_optim_sanity_checks(sdrep, report)
        # })
        convergence_message <- capture.output(
            tryCatch({
                sdrep = TMB::sdreport(em_model_obj)
                SPoRC::post_optim_sanity_checks(sdrep, report)
            }, error = function(e){
                # cat(e$message)
            }),
            type="message"
        )

        convergence_message <- ifelse(length(convergence_message > 1), convergence_message[1], convergence_message)

        # possible error types: Inf in LL, gradient, hessian, non-finite SEs, SEs big, correlations
        if(grepl("Found Inf in joint log-likelihood", convergence_message)){
            parameter_name <- "Unknown"
            problem <- "Inf in JLL"
        }else if(grepl("absolute gradient", convergence_message)){
            split1 <- str_split(convergence_message, "\\. ")[[1]]
            split1 <- split1[-length(split1)]

            s <- str_split(split1, "Parameter: ")[[1]][2]
            s2 <- str_split(s, " ")[[1]]

            parameter_name <- s2[1]
            problem <- paste(s2[3:6], collapse=" ")

        }else if(grepl("Hessian", convergence_message)){
            parameter_name <- "Unknown"
            problem <- "Hessian not positive definite"
        }else if(grepl("non finite elements", convergence_message)){
            parameter_name <- "Unknown"
            problem <- "Non finite SE"
        }else if(grepl("which was greated than tolerance", convergence_message)){

            split1 <- str_split(convergence_message, "\\. ")[[1]]
            split1 <- split1[-length(split1)]

            s <- str_split(split1, "Parameter: ")[[1]][2]
            s2 <- str_split(s, " ")[[1]]

            parameter_name <- s2[1]
            problem <- paste(s2[4:7], collapse = " ")
        }else if(grepl("Parameter pairs:", convergence_message)){

            split1 <- str_split(convergence_message, "\\. ")[[1]]
            split1 <- split1[-length(split1)]

            s <- str_split(split1, "Parameter pairs: ")[[1]][2]
            s2 <- str_split(s, " ")[[1]]

            parameter_name <- paste(s2[1:3], collapse=" ")
            problem <- paste(s2[(length(s2)-2):length(s2)], collapse=" ")
        }

        
        convergence_table$parameter[i] <- parameter_name
        convergence_table$problem[i] <- problem #problem
    }
    return(convergence_table)
}

#' Find Specific MSE EM Model Run from Saved Files
#' 
#' @param hcr_filter name of HCR to find
#' @param om_filter name of OM to find
#' @param simulation_seed simulation seed to find
#' 
#' @return list with mse object for specified model run
#' 
find_model_run <- function(hcr_filter, om_filter, simulation_seed, om_list, seed_list){

    # Identify correct file for OM and seed
    oms <- lapply(om_list, function(x) x$name)
    om_num <- which(om_filter == oms)-1

    max_seed_num <- ceiling(length(seed_list)/22)
    seed_num <- which(simulation_seed == seed_list)
    seed_file <- ceiling(seed_num/22) # 22 seeds per file by default

    file_suffix <- om_num*max_seed_num+seed_file

    # fs <- list.files(file.path(here::here(), "data", "active"), full.names = TRUE)
    # if(!is.null(hcr_filter))
    #     fs <- unlist(sapply(hcr_filter, \(x) fs[grepl(paste0(sub("|", "\\|", sub("/", "", sub(" ", "_", tolower(x))), fixed=TRUE), file_suffix), fs)]))#"_\\d+"), fs)]))
    fs <- c(file.path(here::here(), "data", "active", paste0("mse_runs_", sub("/", "", sub(" ", "_", tolower(hcr_filter))), "_",file_suffix,".RDS")))

    for(f in 1:length(fs)){
        x <- fs[f]
        # print(x)
        m <- readRDS(x)
        mse <- m$mse_objects
        model_run <- list(mse[[length(mse)]])[[1]]
        # rm(mse)
        # print("Loaded model run object")
        if(!(model_run$om$name %in% om_filter)){
            print("Wrong OM, skipping")
            next;
        }
        
        # n_proj_years <- m$mse_options_list$mse_options$n_proj_years
        # n_spinup_years <- m$mse_options_list$mse_options$n_spinup_years
        seeds <- m$seeds
        nsims <- length(seeds)

        if(!(simulation_seed %in% seeds)){
            print("Seed not in this file, skipping")
            next;
        }
        # simulation_number <- which(seeds == simulation_seed)
        print(paste("Found model run in file:", x))
        return(afscOM::listN(m, model_run))
    }
}


compute_dynamic_value <- function(landings, min_price_age, max_price_age, breakpoints=c(15, 30)){
    if(landings < breakpoints[1]){
        return(max_price_age)
    }else if(landings >= breakpoints[1] & landings <= breakpoints[2]){
        return(
            min_price_age + (breakpoints[2]-landings)/(breakpoints[2]-breakpoints[1])*(max_price_age-min_price_age)
        )
    }else{
        return(min_price_age)
    }
}

set_prices <- function(){
    price_age_f_low <- c(0.597895623, 1.320303448, 1.320303448, 1.856562267, 2.610111345, 2.610111345, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 6.01401531, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875, 7.435514875)
    price_age_m_low <- c(0.597895623, 0.597895623, 1.320303448, 1.320303448, 1.856562267, 1.856562267, 1.856562267, 1.856562267, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345, 2.610111345)
    price_data_low <- matrix(c(price_age_f_low, price_age_m_low), nrow=length(price_age_f_low), ncol=2)
    dimnames(price_data_low) <- list("age"=2:31, "sex"=c("F", "M"))

    price_age_f_max <- c(7.917460094, 8.40756497, 8.40756497, 9.944657109, 11.46480347, 11.46480347, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 12.97470389, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658, 14.86275658)
    price_age_m_max <- c(7.917460094, 7.917460094, 8.40756497, 8.40756497, 9.944657109, 9.944657109, 9.944657109, 9.944657109, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347, 11.46480347)
    price_data_max <- matrix(c(price_age_f_max, price_age_m_max), nrow=length(price_age_f_max), ncol=2)
    dimnames(price_data_max) <- list("age"=2:31, "sex"=c("F", "M"))

    return(afscOM::listN(price_data_low, price_data_max))
}
