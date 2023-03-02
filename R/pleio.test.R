#' A function to perform an IV-specific test of nonzero horizontal pleiotropy for many IVs at once
#' 
#' @param Ahat a mxp matrix (p exposures) of standardized GWAS estimates for the m SNPs
#' @param Bhat a mxq matrix (q outcomes) of standardized GWAS estimates for the m SNPs
#' @param Thetahat a pxq matrix of causal effect estimates from MRBEE
#' @param CovThetahat An estimate of the variance-covariance matrix of Thetahat
#' @param SigmaUU an pxpxm array of variance-covariance matrices of measurement errors for the p exposures
#' @param SigmaUV an pxqxm array of covariance matrices of measurement errors for the p exposures and q outcomes
#' @param SigmaVV an qxqxm array of variance-covariance matrices of measurement errors for the q outcomes
#' @param Outliers vector of indices of items to remove from Ahat, Bhat, SigmaUU, SigmaUV, and SigmaVV from their m-index
#' @keywords IMRP pleiotropy test for many SNPs
#' @export 
#' @examples
#' pleio.test()

pleio.test=function(Ahat, Bhat, Thetahat, CovThetahat, SigmaUU, SigmaUV, SigmaVV, Outliers=NULL) {
	.dat=subset.all(Ahat, Bhat, SigmaUU, SigmaUV, SigmaVV, Outliers=Outliers)
	Ahat=as.matrix(.dat$A); Bhat=as.matrix(.dat$B); SigmaUU=.dat$SigmaUU; SigmaUV=.dat$SigmaUV; SigmaVV=.dat$SigmaVV
	m=nrow(Ahat);p=ncol(Ahat);q=ncol(Bhat)
	Si=sapply(1:m, function(h) {
		kronI=kronecker(Ahat[h,], diag(length(Bhat[h,])))
		Upsilonhati=SigmaVV[,,h]+t(Thetahat)%*%SigmaUU[,,h]%*%Thetahat+t(kronI)%*%CovThetahat%*%kronI-2*t(Thetahat)%*%SigmaUV[,,h]
		Si=fS(Bhat[h,], Thetahat, Ahat[h,], Upsilonhati)
	})
	Pvalues=1-pchisq(Si, p*q)
	return(list("Stats"=Si, "Ps"=Pvalues))
}
