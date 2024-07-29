#' Run MSE Simulations Across OMs and HCRs
#' 
#' Wrapper function around `run_mse` and `run_mse_parallel` that
#' allows for simple running of MSE simulations across multiple 
#' OMs and HCRs (and combinations of the two). The wrapper relies
#' on using `expand.grid` to create a factorial combination of OMs
#' and HCRs and runs MSE simulations across all such combinations. 
#'
#' @param om_list named list of operating model list objects
#' @param hcr_list named list of harvest control rule list objects
#' @param seed_list vector of random seeds
#' @param mse_options list of mse_options objects to apply to each
#' OM. Must be either length 1, or length == length(om_list).
#' @param ... additional parameters to pass to the `run_mse` call
#'
#' @return list of MSE simualtion results
#' @export run_mse_multiple
#' 
#' @example
#'
run_mse_multiple <- function(om_list, hcr_list, seed_list, mse_options_list, nyears){
    
    # if(length(mse_options_list) != 1 && length(mse_options_list) != length(om_list)){
    #     stop("Invalid input for parameter `mse_options`. Parameter must have same length as `om_list` or length 1. If length 1, the same set of options will be used across OMs.")
    # }

    # if(length(mse_options_list) == 1){
    #     mse_options_list <- rep(mse_options_list, length(om_list))
    # }

    mse_run_grid <- expand.grid(om=names(om_list), hcr=names(hcr_list), opt=names(mse_options_list))
    mse_objects <- list()

    nsims <- length(seed_list)

    for(i in 1:nrow(mse_run_grid)){
        # omi <- which(names(om_list) == mse_run_grid[i,1])[1]
        om <- om_list[[mse_run_grid[i,1]]]
        hcr <- hcr_list[[mse_run_grid[i,2]]]
        opt <- mse_options_list[[mse_run_grid[i,3]]]

        mse_run <- run_mse_parallel(nsims, seed_list, om, hcr, mse_options=opt, nyears=nyears)
        mse_objects[[i]] <- mse_run

    }

    return(mse_objects)

}
