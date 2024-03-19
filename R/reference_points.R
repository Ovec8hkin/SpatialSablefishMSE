compute_sbpr <- function(nages, mort, mat, waa, sel, ret, F){
    naa <- rep(NA, nages)
    naa[1] <- 1
    zaa <- mort + sel*ret*F
    for(a in 2:(nages-1)){
        naa[a] <- naa[a-1]*exp(-zaa[a-1])
    }
    ssb <- sum(naa*mat*waa, na.rm = TRUE)
    return(ssb)
}

compute_spr <- function(nages, mort, mat, waa, sel, ret, F){
    ssb_unfished <- compute_sbpr(nages, mort, mat, waa, sel, ret, F=0)
    ssb_fished   <- compute_sbpr(nages, mort, mat, waa, sel, ret, F)
    return(ssb_fished/ssb_unfished)
}

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

compute_bx <- function(nages, mort, mat, waa, sel, ret, F, avg_rec){
  return(avg_rec*compute_sbpr(nages, mort, mat, waa, sel, ret, F))
}

calculate_ref_points <- function(nages, mort, mat, waa, sel, ret, avg_rec){
  F35 <- spr_x(nages, mort, mat, waa, sel, ret ,target_x=0.35)
  F40 <- spr_x(nages, mort, mat, waa, sel, ret, target_x=0.40)
  B40 <- compute_bx(nages, mort, mat, waa, sel, ret, F=F40, avg_rec=avg_rec)
  return(list(F40=F40, F35=F35, B40=B40))
}

# assessment <- dget("data/sablefish_assessment_2023.rdat")

# nages <- 30
# mort <- assessment$M[1,1]
# mat <- assessment$growthmat[,"mage.block1"]
# waa <- assessment$growthmat[,"wt.f.block1"]
# sel <- as.vector(assessment$agesel[,"fish5sel.f"])+as.vector(assessment$agesel[,"fish3sel.f"])
# sel <- sel/max(sel)
# ret <- 1
# avg_rec <- mean(assessment$natage.female[,1])


# compute_sbpr(nages, mort, mat, waa, sel, ret, F=0)
# compute_spr(nages, mort, mat, waa, sel, ret, F=0.0957209)
# F40 <- spr_x(nages, mort, mat, waa, sel, ret, target_x=0.40)

# B40 <- compute_bx(nages, mort, mat, waa, sel, ret, F=F40, avg_rec=avg_rec)

# sum(assessment$natage.female[64,]*mat*waa)/B40
