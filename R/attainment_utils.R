#' Stairstep Attainment Function
#' 
#' Defines attainment level based on a three stepped function,
#' with optional phase ins periods where a linear combination of
#' stairsteps is applied. Attainment is assumed to be a function
#' of the established TAC.
#'
#' @param v a total allowable catch (TAC) level
#' @param breakpoints TAC levels that define where a transition
#' in attainment level occurs
#' @param levels attainment levels associated with each set of breakpoints
#' (0 - breakpoint[1], breakpoinut[1]-breakpoint[2], breakpoint[2] - Inf)
#' @param phase_ins optional phase_in width. When specified, a phase ins
#' begin talking effect at phase_ins/2 before each breakpoints value and
#' end at phase_ins/2 after each breakpoint value.
#'
#' @export stairstep_attainment
#'
#' @example \dontrun{
#'    TAC <- 25000
#'    att <- stairstep_attainment(
#'       v = TAC,
#'       breakpoints = c(20000, 30000),
#'       levels = c(0.90, 0.80, 0.65),
#'       phase_ins = 0
#'    )
#' }
#'
stairstep_attainment <- function(v, breakpoints, levels, phase_ins){
   new_breakpoints <- as.vector(sapply(breakpoints, \(x) c(x-phase_ins/2, x+phase_ins/2)))

   if(is.na(v)){
      return(NA)
   }

   if(v < new_breakpoints[1]){
      return(levels[1])
   }else if(v > new_breakpoints[2] & v < new_breakpoints[3]){
      return(levels[2])
   }else if(v > new_breakpoints[4]){
      return(levels[3])
   }else if(v >= new_breakpoints[1] & v <= new_breakpoints[2]){
      return(
         ((levels[2]-levels[1])/(new_breakpoints[2]-new_breakpoints[1]))*(v-new_breakpoints[1]) + levels[1]
      )
   }else if(v >= new_breakpoints[3] & v <= new_breakpoints[4]){
      return(
         ((levels[3]-levels[2])/(new_breakpoints[4]-new_breakpoints[3]))*(v-new_breakpoints[3]) + levels[2]
      )
   }

}

#' Average Region-Fleet TAC Utilization
#' 
#' Average Region-Fleet TAC utilization rates from 2000-2020 as derived
#' from historical regional ABCs, fleet allocations, and region-specific
#' landings by fleet. See `dev/historical_TAC_analyses.R`.
#' 
#' @export average_regflt_tac_utilization
#' 
average_regflt_tac_utilization <- array(
    c(0.644, 0.511, 0.947, 1.03, 1.01, 
      0.185, 0.104, 0.293, 0.835, 0.491), # Fixed and Trawl fleets
    dim = c(5, 2),
    dimnames = list(
        "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"),
        "fleet"=c("Fixed", "Trawl")
    )
)

#' Full TAC Utilization
#' 
#' Fully TAC utilization by the fleets within each region
#' 
#' @export full_utilization
#' 
full_utilization <- array(1, dim=c(5, 2), dimnames=list(
   "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"),
   "fleet"=c("Fixed", "Trawl")
))
