#' Estimate Causal Effect with MRBEE
#'
#' This function estimates the causal effect using a bias-correction estimating equation with an intercept, considering potential pleiotropy and measurement errors.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The last one should be the outcome.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to 30.
#' @param max.eps Tolerance for stopping criteria. Defaults to 1e-4.
#' @param pv.thres P-value threshold in pleiotropy detection. Defaults to 0.05.
#' @param var.est Method for estimating the variance of residual in pleiotropy test. Can be "robust", "variance", or "ordinal". Defaults is robust that estimates the variance of residual using median absolute deviation (MAD).
#' @param FDR Logical. Whether to apply the FDR to convert the p-value to q-value. Defaults to TRUE.
#' @param adjust.method Method for estimating q-value. Defaults to "Sidak".
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to 3.
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy
#' @importFrom MASS rlm
#' @export
#'
MRBEE.IMRP.Egger=function(by,bX,byse,bXse,Rxy,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="robust",FDR=T,adjust.method="Sidak",maxdiff=3){
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
bX=cbind(1,bX)
bXse=cbind(0,bXse)
Rxy=as.matrix(Matrix::bdiag(0,Rxy))
n=length(by)
p=ncol(bX)
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
########## Initial Estimation ############
fit=MASS::rlm(by~bX-1)
theta.ini=fit$coefficient
theta=theta.ini
theta1=10000
e=c(by-bX%*%theta)
indvalid=which(abs(e)<=3*mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
########## Iteration ###################
error=sqrt(sum((theta-theta1)^2))
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-bX%*%theta)
pv=imrpdetect(x=e,theta=theta,RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
Hinv=t(bX[indvalid,])%*%bX[indvalid,]-Rxysum[1:p,1:p]
Hinv=solve(Hinv)
theta=Hinv%*%(t(bX[indvalid,])%*%by[indvalid]-Rxysum[1+p,1:p])
if((norm(theta,"2")/norm(theta.ini,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini,"2")
}
iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}
adjf=n/(length(indvalid)-dim(bX)[2])
D=bX[indvalid,]%*%(Hinv%*%t(bX[indvalid,]))
D=(rep(1,length(indvalid))-diag(D))
D[which(D<0.25)]=0.25
E=-bX[indvalid,]*(e[indvalid]/D)
for(i in 1:length(indvalid)){
E[i,]=E[i,]+RxyList[indvalid[i],p+1,1:p]-c(RxyList[indvalid[i],1:p,1:p]%*%theta)
}
V=t(E)%*%E
covtheta=(Hinv%*%V%*%Hinv)*adjf
r=c((by-bX%*%theta))*byse1
r[indvalid]=0
names(theta)=colnames(bX)
names(r)=rownames(bX)
A=list()
A$theta=theta
A$covtheta=covtheta
A$delta=r
return(A)
}
