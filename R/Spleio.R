#' A wrapping function to perform an IV-specific test of nonzero horizontal pleiotropy for many IVs at once
#' 
#' @param prepDataList the direct output of prepData(), a list of length 5.
#' @param Outliers vector of indices of IVs to remove
#' @export 
#' @examples
#' Spleio()

Spleio=function(prepDataList, Thetahat, CovThetahat, Outliers=NULL) {
  bx=prepDataList$betaX;by=prepDataList$betaY;
  UU=prepDataList$UU;UV=prepDataList$UV;VV=prepDataList$VV
  p=dim(UU)[1];q=dim(VV)[1]
  SigmaUU0=lapply(array.to.list(UU), function(h) rbind(rep(0, p+1), cbind(rep(0, p), h)))
  SigmaUV0=lapply(array.to.list(UV), function(h) rbind(rep(0, q), matrix(h,nr=p,nc=q)))
  SigmaUU0=simplify2array(SigmaUU0); SigmaUV0=simplify2array(SigmaUV0);
  out=pleio.test(cbind(1,bx),by,Thetahat,CovThetahat,SigmaUU0,SigmaUV0,VV)
  return(out)
}
