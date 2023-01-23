#' A function to perform non-iterative MRBEE with added intercept terms
#'
#' @param Ahat mxp matrix of standardized GWAS estimates for p exposures and m IVs
#' @param Bhat mxq matrix of standardized GWAS estimates for q outcomes and m IVs
#' @param SigmaUU pxpxm array of measurement error variance-covariance matrices for the p exposures and m IVs.
#' @param SigmaUV pxqxm array of measurement error covariance matrices for the p exposures, q outcomes, and m IVs.
#' @param SigmaVV qxqxm array of measurement error variance-covariance matrices for the q outcomes and m IVs.
#' @param Outliers vector of indices of items to remove from Ahat, Bhat, SigmaUU, SigmaUV, and SigmaVV from their m-index
#' @keywords MRBEE with intercept
#' @export
#' @examples
#' EE.estimator.intercept()

EE.estimator.intercept=function(Ahat, Bhat, SigmaUU, SigmaUV, SigmaVV, Outliers=NULL) {
  if(length(dim(SigmaVV)) != 3) stop("please enter variance-covariance matrices as 3D arrays")
  .dat=subset.all(Ahat, Bhat, SigmaUU, SigmaUV, SigmaVV, Outliers=Outliers)
  Ahat=as.matrix(.dat$A); Bhat=as.matrix(.dat$B); SigmaUU=.dat$SigmaUU; SigmaUV=.dat$SigmaUV; SigmaVV=.dat$SigmaVV
  m=dim(SigmaVV)[3]; p=ncol(Ahat); q=ncol(Bhat)
  # initial estimate - add intercept terms
  Ahat0=cbind(1, Ahat)
  SigmaUU0=lapply(array.to.list(SigmaUU), function(h) rbind(rep(0, p+1), cbind(rep(0, p), h)))
  SigmaUV0=lapply(array.to.list(SigmaUV), function(h) rbind(rep(0, q), matrix(h,nr=p,nc=q)))
  # back to arrays
  SigmaUU0=simplify2array(SigmaUU0); SigmaUV0=simplify2array(SigmaUV0);
  fit0=EE.estimator(Ahat0, Bhat, SigmaUU0, SigmaUV0, SigmaVV)
  rownames(fit0$Est)=c("intercept", paste0("X", 1:p)); colnames(fit0$Est)=paste0("Y", 1:q)
  # also need to perform pleiotropy test on all SNPs at first iteration
  pleioPs=pleio.test(Ahat0, Bhat, fit0$Est, fit0$Var2, SigmaUU0, SigmaUV0, SigmaVV, Outliers=Outliers)$P
  # Can also perform a global test with this information
  intMat=matrix(1:((p+1)*q),nr=p+1,nc=q); intInds=c(intMat[1,])
  C=matrix(0, nr=q, nc=(p+1)*q); for(i in 1:q) C[i,intInds[i]]=1
  Q=fQ(C, fit0$Est, fit0$Var2); PQ=1-pchisq(Q, q)
  return(list(Est=fit0$Est, Var=fit0$Var2, pleioPs=pleioPs, Q=Q, PQ=PQ,
              SigmaUU0=SigmaUU0,SigmaUV0=SigmaUV0))
}