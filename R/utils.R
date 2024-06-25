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
