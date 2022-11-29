#' A function to perform MR using MRBEE
#'
#' This function performs multivariable/multivariate MR using MRBEE.
#' @param prepDataList the direct output of prepData(), a list of length 5.
#' @param PleioPThreshold The P-value threshold below which a specific IV will be considered horizontally pleiotropic and removed from causal estimation using pleio_test() (see paper).
#' @param FDR logical. If you want to use Benjamini-Hochberg FDR correction instead of defining a P-value threshold for determing horizontal pleiotropy evidence (as in the argument PleioPThreshold), set to TRUE.
#' @param FDR.alpha the uncorrected Type I error rate to be corrected by FDR.
#' @keywords preparing data
#' @export 
#' @examples
#' MRBEE()
 
MRBEE=function(prepDataList, PleioPThreshold=0.05, FDR=FALSE, FDR.alpha=0.05) {
  bx=prepDataList$betaX;by=prepDataList$betaY;
  UU=prepDataList$UU;UV=prepDataList$UV;VV=prepDataList$VV
  it=iter.estimator(A=bx,B=by,SigmaUU=UU,SigmaUV=UV,SigmaVV=VV,PleioPThreshold=PleioPThreshold,FDR=FDR,FDR.alpha=FDR.alpha,Outliers=NULL)
  return(it)
}