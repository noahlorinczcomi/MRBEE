#' A helping function to perform MR using MRBEE
#'
#' MRBEE is intended to be performed using the MRBEE() function, which calls the function below.
#' @param A mxp matrix of standardized exposure GWAS estimates
#' @param B mxq matrix of standardized exposure GWAS estimates
#' @param UU pxpxm array of measurement error variance-covariance matrices for the p exposures and m IVs.
#' @param UV pxqxm array of measurement error covariance matrices for the p exposures, q outcomes, and m IVs.
#' @param VV qxqxm array of measurement error variance-covariance matrices for the q outcomes and m IVs.
#' @keywords iterative MRBEE with IMRP adjustment
#' @export 
#' @examples
#' iter.estimator()

iter.estimator=function(A, B, SigmaUU, SigmaUV, SigmaVV, PleioPThreshold, FDR=FALSE, FDR.alpha=0.05, Outliers=NULL,warn=TRUE) {
  .dat=subset.all(A, B, SigmaUU, SigmaUV, SigmaVV, Outliers=Outliers)
  Ahat=as.matrix(.dat$A); Bhat=as.matrix(.dat$B); SigmaUU=.dat$SigmaUU; SigmaUV=.dat$SigmaUV; SigmaVV=.dat$SigmaVV
  m=nrow(Bhat); p=ncol(Ahat); q=ncol(Bhat)
  fit0=EE.estimator.intercept(Ahat, Bhat, SigmaUU, SigmaUV, SigmaVV, Outliers=Outliers)
  SigmaUU0=fit0$SigmaUU0;SigmaUV0=fit0$SigmaUV0;
  initial_ests=fit0$Est; initial_vcov=fit0$Var
  pleioPs0=fit0$pleioPs
  Q0=fit0$Q; PQ0=1-pchisq(Q0,q)
  intMat=matrix(1:((p+1)*q),nr=p+1,nc=q); intInds=c(intMat[1,])
  C=matrix(0, nr=q, nc=(p+1)*q); for(i in 1:q) C[i,intInds[i]]=1  

  k=0; thetadiff=1; tdd=1
  if(FDR) Outliers=suppressWarnings(which(stats::p.adjust(pleioPs0,method="BH")<FDR.alpha)) else Outliers=which(pleioPs0<PleioPThreshold)
  diffs=numeric(); PQiter=PQ0
  while(k<50 & thetadiff>(0.0001*(p*q+p)) & tail(tdd,1)!=0 & length(Outliers)<dim(SigmaUU)[3]) { # times that term because I want it to be sensitive to many phenotypes
    k=k+1
    fitk=EE.estimator.intercept(Ahat, Bhat, SigmaUU, SigmaUV, SigmaVV, Outliers)
    pleioPs0=pleio.test(cbind(1,Ahat), Bhat, as.matrix(fitk$Est), fitk$Var, SigmaUU0, SigmaUV0, SigmaVV, Outliers=NULL)$Ps
    if(FDR) Outliers=suppressWarnings(which(p.fdr(pleioPs0,just.fdr=TRUE)<FDR.alpha)) else Outliers=which(pleioPs0<PleioPThreshold)
    thetadiff=sum(abs(fitk$Est-fit0$Est))
    fit0=fitk
    fitk_int=EE.estimator.intercept(Ahat, Bhat, SigmaUU, SigmaUV, SigmaVV, Outliers)
    Qk=fQ(C, fitk_int$Est, fitk_int$Var)
    PQk=1-pchisq(Qk, q)
    diffs[k]=thetadiff
    if(k>1) tdd=c(tdd, diffs[k]-diffs[k-1]) # just making sure we aren't iterating the same results over and over
  }
  if(warn & PQk<0.05) warning("the final global test of unbalanced pleiotropy has P<0.05. If GWAS sample sizes are small (e.g. < ~50k), consider setting PleioPThreshold to a larger value (e.g., 0.05)")
  # done
  out=list(
    CausalEstimates=fitk$Est,
    VCovCausalEstimates=fitk$Var,
    CausalEstimatePS=2*pnorm(-abs(fitk$Est/sqrt(diag(fitk$Var)))),
    InitialCausalEstimates=initial_ests,
    VcovInitialCausalEstimates=initial_vcov,
    PleiotropyIndices=Outliers,
    nIterations=k,
    Qstart=Q0,
    Qend=Qk,
    ThetaDifferences=diffs
  )
  return(out)
}
