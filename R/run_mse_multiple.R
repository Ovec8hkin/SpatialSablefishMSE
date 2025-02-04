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
run_mse_multiple <- function(om_list, hcr_list, seed_list, mse_options_list, nyears, diagnostics=FALSE, save=FALSE){
    
    # if(length(mse_options_list) != 1 && length(mse_options_list) != length(om_list)){
    #     stop("Invalid input for parameter `mse_options`. Parameter must have same length as `om_list` or length 1. If length 1, the same set of options will be used across OMs.")
    # }

    # if(length(mse_options_list) == 1){
    #     mse_options_list <- rep(mse_options_list, length(om_list))
    # }

    mse_run_grid <- expand.grid(om=names(om_list), hcr=names(hcr_list), opt=names(mse_options_list))
    mse_objects <- list()

    nsims <- length(seed_list)

    max_sims <- 66
    nsim_iters <- ifelse(nsims < max_sims, 1, round(nsims / max_sims)+1)

    m <- 1
    counter <- 0
    for(i in 1:nrow(mse_run_grid)){
        # omi <- which(names(om_list) == mse_run_grid[i,1])[1]
        om <- om_list[[mse_run_grid[i,1]]]
        hcr <- hcr_list[[mse_run_grid[i,2]]]
        opt <- mse_options_list[[mse_run_grid[i,3]]]
    
        for(j in 1:nsim_iters){
            seeds <- seed_list[(((j-1)*max_sims)+1):min((j*max_sims), length(seed_list))]
            nsims2 <- length(seeds)
            mse_run <- run_mse_parallel(nsims2, seeds, om, hcr, mse_options=opt, nyears=nyears, diagnostics=diagnostics)
            
            counter <- counter+1
            mse_objects[[j]] <- mse_run

            print(save)
            if(save || nsim_iters>1){
                # Going to save files by HCR
                filename <- file.path(here::here(), "data", "active", paste0("mse_runs_", sub("/", "", sub(" ", "_", tolower(hcr$name))), "_",counter,".RDS"))
                obj <- listN(mse_objects, om, hcr, mse_options_list, seeds)
                saveRDS(obj, file=filename)
                # m <- m+1
                mse_objects <- list()
            }
        }
        
        next_hcr <- hcr_list[[mse_run_grid[i+1, 2]]]
        if(!identical(next_hcr, hcr)){
            counter <- 0
        }

    }

    return(mse_objects)

}