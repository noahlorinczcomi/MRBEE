#' A function to prepare GWAS data for MRBEE
#'
#' This function prepares GWAS data for MR with MRBEE.
#' @param betaX a mxp matrix of standardized GWAS estimates for p exposures and m genetic instrumental variables.
#' @param betaY a mxq matrix of standardized GWAS estimates for q outcomes and m genetic instrumental variables.
#' @param seX a mxp matrix of standardized GWAS standard errors for p exposures and m genetic instrumental variables (columns match order of betaX columns).
#' @param seY a mxp matrix of standardized GWAS standard errors for p outcomes and m genetic instrumental variables (columns match order of betaY columns).
#' @param R (p+q)x(p+q) matrix of correlations between GWAS estimates. See corrMatrix.py software.
#' @keywords preparing data
#' @export
#' @examples
#' prepData()

prepData=function(betaX,betaY,seX,seY,R,verbose=TRUE) {
  if(verbose) cat("NOTE: this software assumes the following:\n(i) the ith row of all objects passed to it corresponds to the ith IV\n(ii) all values passed are in the standardized scale of your choice,\n(iii) all IVs are independent\n")
  betaX=as.matrix(betaX);betaY=as.matrix(betaY)
  seX=as.matrix(seX);seY=as.matrix(seY);R=as.matrix(R)
  if(any(dim(betaX)!=dim(seX)) | any(dim(betaY)!=dim(seY))) stop("dimensions of betas and standard errors do not all match")
  if(!all(colnames(betaX)==colnames(seX))) stop("colnames of betaX do not match colnames of seX. Colnames must match and have the same order.")
  if(!all(colnames(betaY)==colnames(seY))) stop("colnames of betaY do not match colnames of seY. Colnames must match and have the same order.")
  if(dim(R)[1]!=(ncol(betaX)+ncol(betaY))) stop("dimensions of R do not match number of exposures and outcomes provided")
  cnR=colnames(R);cnX=colnames(betaX);cnY=colnames(betaY);cnXs=colnames(seX);cnYs=colnames(seY)
  if(!all(cnR %in% rownames(R))) stop("row and column names of R do not match each other")
  if(!(all(c(cnX,cnY,cnXs,cnYs) %in% cnR))) stop("please make column names and order for betaX match seX and for betaY match seY")
  indY=ff(R,cnY); indX=ff(R,cnX); R=R[c(indY,indX),c(indY,indX)]
  m=nrow(betaX);p=ncol(betaX);q=ncol(betaY)
  Ryy=R[1:q,1:q]; Rxx=R[(q+1):(p+q),(q+1):(p+q)];Rxy=R[(q+1):(q+p),1:q]
  UU=array(dim=c(p,p,m));UV=array(dim=c(p,q,m));VV=array(dim=c(q,q,m))
  for(i in 1:m) {
    sx=seX[i,];sy=as.matrix(seY[i,])
    UU[,,i]=diag(sx)%*%Rxx%*%diag(sx)
    VV[,,i]=diag(sy)%*%Ryy%*%diag(sy)
    UV[,,i]=diag(sx)%*%Rxy%*%diag(sy)
  }
  prepDataList=list(betaX=betaX,betaY=betaY,UU=UU,UV=UV,VV=VV)
  return(prepDataList)
}
