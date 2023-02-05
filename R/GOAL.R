#'@title Genelized outcome-adaptive LASSO estimation
#'
#' @param data Data that include outcome, treatment and all of covariates
#' @param var.list names of potential unknown confounders
#' @param covar names of known confounders
#' @param Trt names of treatment
#' @param out names of outcome
#' @param lambda_vec  a vector of possible lambda values and defualt is c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
#' @param gamma_convergence the value of gamma converagence and defualt is 2
#'
#' @return ATE is the linear dose-response function estimated using IPTW method.    selectedVar is covariates that GOAL selected. Ideally, it include confounders and prognostic covariates.    lambda is the optimal lambda corresponding to minimized DWC.
#' @export
#'
#' @examples
GOAL<-function (data,var.list,covar=NULL,Trt="Trt",out="Y",lambda_vec=c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),gamma_convergence=2) {
  Data<-data
  n<-dim(Data)[1]
  #??׼??
  temp.mean <- colMeans(Data[,var.list])
  Temp.mean <- matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] <- Data[,var.list] - Temp.mean
  temp.sd <- apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd <- matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] <- Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  #outcome model
  if (is.null(covar)) {
    allvar<-var.list
  } else {
    allvar<-c(covar,var.list)
  }

  y.form <- formula(paste(out,"~",paste(c(Trt,allvar),collapse="+")))
  lm.Y <- lm(y.form,data=Data)
  betaXY <- coef(lm.Y)[var.list]
  lambda_vec <- lambda_vec
  names(lambda_vec) <- as.character(lambda_vec)
  gamma_convergence_factor <- gamma_convergence
  gamma_vals <- 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) <- names(lambda_vec)

  wAMD_vec=rep(NA, length(lambda_vec))
  ATE_lambda=vector(mode = "list", length(lambda_vec))
  names(ATE_lambda)=names(wAMD_vec)=names(lambda_vec)
  coeff_XA <- as.data.frame(matrix(NA,nrow=length(var.list),
                                   ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) <- var.list
  CBPS.var<-list()
  w.full.form <- formula(paste(Trt,"~",paste(allvar,collapse="+")))
  for(lil in names(lambda_vec)){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    if (is.null(covar)) {
      oal_pen <- adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig))
      logit_oal <- lqa.formula(w.full.form, data=Data, penalty=oal_pen,
                               family=gaussian())
      coeff_XA[var.list,lil] <- coef(logit_oal)[var.list]

      CBPS.var[[lil]]<-c(rownames(coeff_XA[var.list,])[which(round(coeff_XA[var.list,lil],5)!=0)])
      if (length(CBPS.var[[lil]])!=0) {
        w.model<-formula(paste(Trt,"~",paste(CBPS.var[[lil]],collapse="+")))
        CBPS_fit<-npCBPS(w.model,data=Data,corprior=.1/n, print.level=1)
        Data[,paste("w",lil,sep="")]<-CBPS_fit$weights
      } else {
        Data[,paste("w",lil,sep="")]<-1
      }
      wAMD_vec[lil] <- wAMD_function(DataM=Data,varlist=var.list,trt.var=Trt,
                                     wgt=paste("w",lil,sep=""),beta=coef(lm.Y)[var.list])$wAMD
      ATE_lambda[[lil]] <- ATE_est(fY=Data[,out],fw=Data[,paste("w",lil,sep="")],
                                   fA=Data[,Trt])

    } else{
      oal_pen <- adaptive.lasso(lambda=n^(il),al.weights = c(rep(0,length(covar)),abs(betaXY)^(-ig)))
      logit_oal <- lqa.formula(w.full.form, data=Data, penalty=oal_pen,
                               family=gaussian())
      coeff_XA[var.list,lil] <- coef(logit_oal)[var.list]

      CBPS.var[[lil]]<-c(covar,rownames(coeff_XA[var.list,])[which(round(coeff_XA[var.list,lil],5)!=0)])
      if (length(CBPS.var[[lil]])!=0) {
        w.model<-formula(paste(Trt,"~",paste(CBPS.var[[lil]],collapse="+")))
        CBPS_fit<-npCBPS(w.model,data=Data,corprior=.1/n, print.level=1)
        Data[,paste("w",lil,sep="")]<-CBPS_fit$weights
      } else {
        Data[,paste("w",lil,sep="")]<-1
      }
      wAMD_vec[lil] <- wAMD_function(DataM=Data,varlist=c(covar,var.list),trt.var=Trt,
                                     wgt=paste("w",lil,sep=""),beta=coef(lm.Y)[c(covar,var.list)])$wAMD
      ATE_lambda[[lil]] <- ATE_est(fY=Data[,out],fw=Data[,paste("w",lil,sep="")],
                                   fA=Data[,Trt])
    }
  }
  ATE<-ATE_lambda[[which.min(wAMD_vec)]]
  Svar<-CBPS.var[[which.min(wAMD_vec)]]
  lambda<-names(wAMD_vec)[which.min(wAMD_vec)]
  GOAL_results<-list(
    ATE=ATE,
    selectedVar=Svar,
    lambda=lambda
  )
  return (GOAL_results)
}
