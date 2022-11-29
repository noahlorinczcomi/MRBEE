#' A function to perform non-iterative MRBEE
#'
#' @param A mxp matrix of standardized GWAS estimates for p exposures and m IVs
#' @param B mxq matrix of standardized GWAS estimates for q outcomes and m IVs
#' @param SigmaUU pxpxm array of measurement error variance-covariance matrices for the p exposures and m IVs.
#' @param SigmaUV pxqxm array of measurement error covariance matrices for the p exposures, q outcomes, and m IVs.
#' @param SigmaVV qxqxm array of measurement error variance-covariance matrices for the q outcomes and m IVs.
#' @param Outliers vector of indices of items to remove from Ahat, Bhat, SigmaUU, SigmaUV, and SigmaVV from their m-index
#' @param perform.pleio logical indicating whether IV-specific horizontal pleiotropy tests should be performed using the statistic Spleio
#' @export
#' @examples
#' EE.estimator.intercept()
 
EE.estimator=function(A, B, SigmaUU, SigmaUV, SigmaVV, Outliers=NULL, perform.pleio=TRUE) {
  if(!is.array(SigmaUU) | !is.array(SigmaUV) | !is.array(SigmaVV)) stop("please enter variance-covariance matrices as arrays")
  A=as.matrix(A); B=as.matrix(B)
  .dat=subset.all(A, B, SigmaUU, SigmaUV, SigmaVV, Outliers=Outliers)
  A=as.matrix(.dat$A); B=as.matrix(.dat$B); SigmaUU=.dat$SigmaUU; SigmaUV=.dat$SigmaUV; SigmaVV=.dat$SigmaVV
  m=nrow(A); p=ncol(A); q=ncol(B)
  Theta=solve(t(A) %*% A - ssum(SigmaUU)) %*% (t(A) %*% B - ssum(SigmaUV))
  I=diag(q)
  A1theta=1/m*kronecker(I, t(A) %*% A - ssum(SigmaUU))
  B1theta=B2theta=matrix(0,nr=p*q,nc=p*q)
  for(i in 1:m) {
    res=B[i,]-t(Theta) %*% A[i,]
    B1theta=B1theta+kronecker(I, A[i,]) %*% res %*% t(res) %*% kronecker(I, t(A[i,]))
    g=kronecker(I, A[i,]) %*% (B[i,] - kronecker(I, t(A[i,])) %*% c(Theta)) - c(SigmaUV[,,i]) + kronecker(I, SigmaUU[,,i]) %*% c(Theta) # nolint
    B2theta=B2theta + g %*% t(g)
  }
  .Var1=1/m*solve(A1theta) %*% (B1theta/m) %*% t(solve(A1theta))
  .Var2=1/m*solve(A1theta) %*% (B2theta/m) %*% t(solve(A1theta))
  if(perform.pleio) {
    pleioPs=pleio.test(A, B, Theta, .Var2, SigmaUU, SigmaUV, SigmaVV, Outliers=NULL)$Ps
  } else {
    pleioPs=rep(NA,m)
  }
  out=list(Est=Theta, Var=.Var1, Var2=.Var2, pleioPs=pleioPs)
  return(out)
}