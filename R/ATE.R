#' causal dose-response curve
#'
#' @param fY outcome
#' @param fw balance weights
#' @param fA a continuous treatment
#'
#' @return tt_ATE1 is a linear dose-response function estimated using IPTW method.
#' @export
#'
#' @examples
ATE_est = function(fY,fw,fA){
  dataset<- data.frame(fY,fA,fw)
  design.b<-svydesign(ids= ~1, weights = ~fw,data= dataset)
  tt_ATE <-svyglm(fY~fA,design = design.b)
  tt_ATE1<-list(Trt_coef=summary(tt_ATE)$coefficients,ATE_wreg=tt_ATE)
  return(tt_ATE1)
}
