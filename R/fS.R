#' A helping function to perform an IV-specific test of nonzero horizontal pleiotropy
#' 
#' This function is called from within the MRBEE() function
#' @param Bhati a qx1 vector (q outcomes) of standardized GWAS estimates for the ith SNP
#' @param Thetahat a pxq (q outcomes) matrix of causal effect estimates
#' @param Ahati a px1 vector (p exposures) of standardized GWAS estimates for the ith SNP
#' @param Upsilonhati an estimated variance-covariance calculated from within MRBEE(). Any consistent estimator of the variance of Bhati-t(Thetahat)%*%Ahati
#' @keywords IMRP horizontal pleiotropy test
#' @export 
#' @examples
#' fS()

fS=function(Bhati, Thetahat, Ahati, Upsilonhati) t(Bhati-t(Thetahat) %*% Ahati) %*% solve(Upsilonhati) %*% (Bhati-t(Thetahat) %*% Ahati)