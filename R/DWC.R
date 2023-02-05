#' dual-weigt correlation (DWC) for selecting optimal lambda
#'
#' @param DataM Data that include outcome, treatment and all of covariates
#' @param varlist names of potential confounders including known and unknown
#' @param trt.var names of treatment
#' @param wgt balancing weights
#' @param beta the coefficient of covariates in the unpenalized "full" outcome model
#'
#' @return wAMD is a vector of dual-weight correlation.
#' @export
#'
#' @examples
wAMD_function <- function(DataM,varlist,trt.var,wgt,beta){
  diff_vec <- rep(NA,length(beta))
  names(diff_vec) <- varlist
  for(jj in 1:length(varlist)){
    diff_vec[jj]<-abs(weightedCorr(DataM[,trt.var],DataM[,varlist[jj]],
                                   method="Pearson",weights=DataM[,wgt]))
  }
  wdiff_vec = diff_vec * abs(beta)
  wAMD = c(sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret)
}
