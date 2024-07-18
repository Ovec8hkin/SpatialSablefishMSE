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
