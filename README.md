#GOAL: Generalized Outcome Adaptive LASSO  
Implements the generalized outcome-adaptive LASSO (GOAL) proposed by Qian Gao et al. (2021). The GOAL was designed for the unbiased estimation of the dose-response function (DRF) from high-dimensional observation data. Confounders and predictors of outcome are selected by minimizing dual-weight coefficient (DWC), which considered covariates balance. Simultaneously, the balace weights are estimated and applied to estimate DRF using Inverse Probability Weighting (IPW) method.  
##Installation  
library(devtools)  
install_local(“/path/to/lqa_1.0-3.tar.gz”)
install_github("qg2023/GOAL")
##Example
###Simulation Data
n<-200
p<-100
mean_x <- 0 
sig_x <- 1
rho <- 0
Beta<-0
var.list<-paste("X",1:p,sep="")
Sigma <- matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list)) 
diag(Sigma) <- sig_x^2
Mean_x <- rep(mean_x,p)
set.seed(123)
Data <- as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma,empirical = FALSE))
colnames(Data)<-var.list
betaOut<-c(1,1,1,1,0,0)
betaTrt<-c(1,1,0,0,1,1)
TrueVar<-as.matrix(Data[,c(1:6)])
Data$Trt<-TrueVar%*%betaTrt+rnorm(n=n,mean=0,sd=1)
Data$Y<-Beta*Data$Trt+TrueVar%*%betaOut+rnorm(n=n,mean=0,sd=1) 
###GOAL
library(lqa)
library(CBPS)
library(wCorr)
library(survey)
GOAL(data=Data,var.list=var.list,Trt="Trt",out="Y")
