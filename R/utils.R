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
