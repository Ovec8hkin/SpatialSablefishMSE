compute_naapr <- function(nages, mort, mat, waa, sel, ret, F){
    naa <- rep(NA, nages)
    naa[1] <- 1
    zaa <- mort + sel*ret*F
    for(a in 2:(nages)){
        naa[a] <- naa[a-1]*exp(-zaa[a-1])
    }
    return(naa)
}

#' Compute Spawning Biomass per Recruit (SBPR)
#'
#' Compute SBPR under a given level of fishing mortality.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param F instantenous fishing mortality rate
#'
#' @export
#'
#' @example
#'
compute_sbpr <- function(nages, mort, mat, waa, sel, ret, F){
    naa <- compute_naapr(nages, mort, mat, waa, sel, ret, F)
    ssb <- sum(naa*mat*waa, na.rm = TRUE)
    return(ssb)
}

#' Compute Spawning Potential Ratio (SPR)
#' 
#' Compute SPR, the ratio of spawning biomass under a 
#' given leve of fishing mortality relative to unfished
#' spawning biomass.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param F instantenous fishing mortality rate
#'
#' @export
#'
#' @example
#'
compute_spr <- function(nages, mort, mat, waa, sel, ret, F){
    ssb_unfished <- compute_sbpr(nages, mort, mat, waa, sel, ret, F=0)
    ssb_fished   <- compute_sbpr(nages, mort, mat, waa, sel, ret, F)
    return(ssb_fished/ssb_unfished)
}

#' Find F that yields a given SPR%
#' 
#' Use bisection algorithm to identify the level of
#' fishing mortality required to yield an SPR of x%.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param target_x desired SPR proportion
#'
#' @export
#'
#' @example
#'
spr_x <- function(nages, mort, mat, waa, sel, ret, target_x=0.35){
    range <- vector(length=2)
    range[1] <- 0
    range[2] <- 2
    n.iter <- 20
    i <- 1
    for(i in 1:n.iter) {
      midpoint <- mean(range)
      spr <- compute_spr(nages, mort, mat, waa, sel, ret, F=midpoint)
      if(spr > target_x) {
        range[1] <- midpoint
        range[2] <- range[2]
      }else {
        range[1] <- range[1]
        range[2] <- midpoint
      }
    }
    Fx <- midpoint
    return(Fx)
}

#' Compute average SSB under a level of F
#' 
#' Compute average SSB under a given level of
#' fishing mortality.
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param F instantenous fishing mortality rate
#' @param avg_rec average recruitment
#'
#' @export
#'
#' @example
#'
compute_bx <- function(nages, mort, mat, waa, sel, ret, F, avg_rec){
  return(avg_rec*compute_sbpr(nages, mort, mat, waa, sel, ret, F))
}

#' Calculate NPFMC groundfish reference points
#' 
#' Calculate F_OFL (F_35%), F_ABC (F_40%), and B40
#' (SSB under F_ABC).
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param avg_rec average recruitment
#'
#' @return list of F40, F35, and B40
#' @export calculate_npfmc_ref_points
#'
#' @example
#'
calculate_npfmc_ref_points <- function(nages, mort, mat, waa, sel, ret, avg_rec){
  F35 <- spr_x(nages, mort, mat, waa, sel, ret, target_x=0.35)
  F40 <- spr_x(nages, mort, mat, waa, sel, ret, target_x=0.40)
  B40 <- compute_bx(nages, mort, mat, waa, sel, ret, F=F40, avg_rec=avg_rec)
  return(list(F40=F40, F35=F35, B40=B40))
}

#' Calculate SPR Based Reference Points
#'
#' Generalized function for computing and F and biomass based reference point
#' for a given SPR level
#'
#' @param nages number of ages in age structure
#' @param mort instantaneous natural mortality rate
#' @param mat maturity-at-age vector
#' @param waa weight-at-age vector
#' @param sel total selectivity-at-age vector
#' @param ret total retention-at-age vector
#' @param avg_rec average recruitment
#' @param spr_target target SPR level. If length(spr_target) == 1, same
#' target SPR level assumed for the F and B reference points. If 
#' length(spr_target) == 2, element 1 corresponds to the spr target for 
#' F, while element 2 corresponds to the spr target for B.
#'
#' @return list of F reference point and biomass reference point
#' @export calculate_ref_points
#'
#' @example
#'
calculate_ref_points <- function(nages, mort, mat, waa, sel, ret, avg_rec, spr_target){
  if(length(spr_target) == 1){
    spr_target = rep(spr_target, 2)
  }

  Fmax <- spr_x(nages, mort, mat, waa, sel, ret, target_x=0.35) # this is F_OFL which can't be eclipsed
  Fref <- spr_x(nages, mort, mat, waa, sel, ret, target_x=spr_target[1])

  f2 <- spr_x(nages, mort, mat, waa, sel, ret, target_x=spr_target[2])
  Bref <- compute_bx(nages, mort, mat, waa, sel, ret, F=f2, avg_rec=avg_rec)
  B0 <- compute_bx(nages, mort, mat, waa, sel, ret, F=0, avg_rec=avg_rec)
  return(list(Fmax=Fmax, Fref=Fref, Bref=Bref, B0=B0))
}
