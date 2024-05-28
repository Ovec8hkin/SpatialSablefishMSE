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
#' @param ... additional parameters to pass to the `run_mse` call
#'
#' @return list of MSE simualtion results
#' @export run_mse_multiple
#' 
#' @example
#'
run_mse_multiple <- function(om_list, hcr_list, seed_list, ...){
    
    mse_run_grid <- expand.grid(om=om_list, hcr=hcr_list)
    mse_objects <- list()

    nsims <- length(seed_list)

    for(i in 1:nrow(mse_run_grid)){
        om <- om_list[[mse_run_grid[i,1]]]
        hcr <- hcr_list[[mse_run_grid[i,2]]]

        mse_run <- run_mse_parallel(nsims, seed_list, om, hcr, ...)
        mse_objects[[i]] <- mse_run

    }

    return(mse_objects)

}
