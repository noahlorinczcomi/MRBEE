#' A function to calculate the bias-correction terms used in MRBEE
#'
#' @param merged_data GWAS data containing estimates, standard errors, and effect allele columns for all exposures and outcomes (not just restricted to IVs; should include all SNPs avaialable)
#' @param Ests vector of column names in `merged_data` of association estimates (betas) for all outcomes and exposures (outcome names first)
#' @param SEs vector of column names in `merged_data` of standard errors (SEs) for all outcomes and exposures (outcome names first) ordered the same as `Ests`
#' @param effect_alleles vector of column names in `merged_data` of effect alleles for all outcomes and exposures (outcome names first) ordered the same as `Ests` and `SEs`
#' @param harmonise_to index of `Ests`, `SEs`, and `effect_alleles` to harmonise alleles to (default is the phenotype corresponding to the first index in `Ests`,`SEs`,`effect_alleles`, which will harmonise all exposure effect alleles to the outcome effect allele for a single outcome)
#' @param pval_threshold only SNPS with P>this threshold will be used to estimate the correlation between measurement errors for all phenotypes
#' @param verbose logical do you want to be reminded that the ordering of column names in `Ests`, `SEs`, and `effect_alleles` should each correspond to the same ordered list of phenotypes used in MR? (TRUE=yes, FALSE=no)
#' @keywords MRBEE
#' @export
#' @examples
#' biasTerms()

biasTerms=function(merged_data,Ests,SEs,effect_alleles,harmonise_to=1,pval_threshold=0.05,verbose=TRUE) {
  # merged_data: merged outcome and exposure GWAS
  # test_statistics: test statistics for outcome and exposures
  # effect_alleles: names of effect allele columns for outcome and exposures (ordered same as test_statistics)
  # harmonise_to: index in test_statistics,effect_alleles, to harmonise alleles to
  # pval_threshold: only SNPs with P>this threshold will be used to calculate R 
  if(verbose) cat(
    '    order of column names for effect_alleles should be 
    in same order as that for Ests and SEs')
  merged_data=as.data.frame(merged_data) # not data.table
  e=merged_data[,Ests]
  s=merged_data[,SEs]
  acut=merged_data[,effect_alleles]
  toHarm=acut!=acut[,harmonise_to]
  e[toHarm]=-e[toHarm]
  dcut=e/s #
  q=qnorm(1-pval_threshold/2)
  keep=apply(dcut,1,function(h) all(abs(h)<q))
  R=cor(na.omit(dcut[keep,])); colnames(R)=rownames(R)=Ests
  Ncor=sum(keep,na.rm=T)
  list(R=R,Ncor=Ncor, EstHarm=as.matrix(e), SEHarm=as.matrix(s))
}
