# The InitialHaplo class contains raw data from MaCS or another program

#InitialHaplo----
#' @title Initial Haplotypes
#' 
#' @description Contains simulated initial haplotypes
#' 
#' @slot nChr number of chromosomes
#' @slot nHaplo number of haplotypes per subpopulation
#' @slot chrData list containing chromosome haplotypes and genetic maps
#' 
#' @export
setClass("InitialHaplo",
         slots=c(nChr="integer",
                 nHaplo="integer",
                 chrData="list"))
