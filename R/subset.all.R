#' A helping function to subset the set of IVs
#'
#' @param A mxp matrix of standardized GWAS estimates for p exposures and m IVs
#' @param B mxq matrix of standardized GWAS estimates for q outcomes and m IVs
#' @param SigmaUU pxpxm array of measurement error variance-covariance matrices for the p exposures and m IVs.
#' @param SigmaUV pxqxm array of measurement error covariance matrices for the p exposures, q outcomes, and m IVs.
#' @param SigmaVV qxqxm array of measurement error variance-covariance matrices for the q outcomes and m IVs.
#' @param Outliers vector of indices of items to remove from Ahat, Bhat, SigmaUU, SigmaUV, and SigmaVV from their m-index
#' @keywords preparing data
#' @export 
#' @examples
#' subset.all()

subset.all=function(A, B, SigmaUU, SigmaUV, SigmaVV, Outliers) {
	if(is.null(Outliers) | length(Outliers) == 0) {
		A=as.matrix(A); B=as.matrix(B)
		return(list(A=A, B=B, SigmaUU=SigmaUU, SigmaUV=SigmaUV, SigmaVV=SigmaVV))
	} else {
		A=as.matrix(A); A=A[-Outliers,]; B=as.matrix(B); B=B[-Outliers,]
		SigmaUU=SigmaUU[,,-Outliers, drop=FALSE]
		SigmaVV=SigmaVV[,,-Outliers,drop=FALSE]
		SigmaUV=SigmaUV[,,-Outliers,drop=FALSE]
		return(list(A=A, B=B, SigmaUU=SigmaUU, SigmaUV=SigmaUV, SigmaVV=SigmaVV))
	}
}