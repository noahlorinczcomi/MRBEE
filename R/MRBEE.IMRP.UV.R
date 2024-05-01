#' Univariable MR Estimation using MRBEE Method
#'
#' This function performs univariable Mendelian Randomization (MR) analysis using
#' the MRBEE method.
#'
#' @param by Vector of GWAS effect sizes for the outcome (n x 1).
#' @param bx Vector of GWAS effect sizes for the exposure (n x 1).
#' @param byse Vector of standard errors (SE) for the GWAS effect sizes of the outcome (n x 1).
#' @param bxse Vector of SE for the GWAS effect sizes of the exposure (n x 1).
#' @param Rxy Correlation matrix (p+1 x p+1) of the exposures and outcome, with the outcome being the last.
#' @param max.iter Maximum number of iterations for the estimation process. Defaults to 30.
#' @param max.eps Tolerance level for convergence. Defaults to 1e-4.
#' @param pv.thres P-value threshold for pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the standard error in the pleiotropy test. Can be "robust", "variance", or "ordinal".
#' @param FDR Logical indicating whether to apply False Discovery Rate (FDR) correction. Defaults to TRUE.
#' @param adjust.method Method for estimating q-values, defaults to "Sidak".
#' @param var.method Method for estimating variance of causal effect, defaults to "sandwich".
#  @param sampling.time Bootstrap time, defaults to 100.
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy
#' @importFrom MASS rlm
#' @export


MRBEE.IMRP.UV=function(by,bx,byse,bxse,Rxy,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="robust",FDR=T,adjust.method="Sidak",var.method="sandwich",sampling.time=100){
by=by/byse
byseinv=1/byse
bx=bx*byseinv
bxse=bxse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
RxyList=IVweight(byse,bxse,Rxy)
########## Initial Estimation ############
fit=MASS::rlm(by~bx-1)
theta=fit$coefficient
theta1=10000
e=c(by-bx*theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
########## Iteration ###################
error=abs(theta-theta1)
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-bx*theta)
pv=imrpdetect(x=e,theta=theta,RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
h=sum(bx[indvalid]^2)-sum(bxse[indvalid]^2*Rxy[1,1])
g=sum(bx[indvalid]*by[indvalid])-Rxy[1,2]*sum(bxse[indvalid]*byse[indvalid])
theta=g/h
iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}
if(var.method=="sandwich"){
adjf=n/(length(indvalid)-1)
Hat=outer(bx[indvalid],bx[indvalid])/h
Hat=1-diag(Hat)
Hat[Hat<0.5]=0.5
e[indvalid]=e[indvalid]/Hat
E=-bx[indvalid]*e[indvalid]+bxse[indvalid]*byse[indvalid]-bxse[indvalid]^2*theta
vartheta=sum(E^2)/h^2*adjf
}else{
thetavec=c(1:sampling.time)
for(i in 1:sampling.time){
indvalidj=sample(indvalid,length(indvalid),replace=T)
h=sum(bx[indvalidj]^2)-sum(bxse[indvalidj]^2*Rxy[1,1])
g=sum(bx[indvalidj]*by[indvalidj])-Rxy[1,2]*sum(bxse[indvalidj]*byse[indvalidj])
thetavec[i]=g/h
}
vartheta=var(thetavec)*n/(length(indvalid)-1)
}
A=list()
A$theta=theta
A$vartheta=vartheta
r=c(by-bx*theta)*byse1
r[indvalid]=0
names(r)=rownames(bx)
A$delta=r
if(var.method!="sandwich"){
A$sampling.theta=mean(thetavec)
}
return(A)
}
