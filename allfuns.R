#' Calculate bias-correction terms from full GWAS data
#'
#' @param merged_data dataframe of merged exposure(s) and outcome GWAS
#' @param dNames length 3 list of with elements that are vectors of column names in merged_data corresponding to effect sizes, SEs, and effect alleles for all exposures and outcome. Order of elements in list vectors must correspond to the same phenotypes
#' @param harmonise_to phenotype to which effect alleles should be harmonised to
#' @param pval_threshold only SNPs with all P>pval_threshold will be used to calculate bias-correction terms
#' @param verbose TRUE/FALSE should some warnings be printed
#' @keywords 
#' @export
#' @examples
#' biasTerms()

biasTerms=function(merged_data,dNames,harmonise_to=1,pval_threshold=0.05,verbose=TRUE) {
  # merged_data: merged outcome and exposure GWAS
  # test_statistics: test statistics for outcome and exposures
  # effect_alleles: names of effect allele columns for outcome and exposures (ordered same as test_statistics)
  # harmonise_to: index in test_statistics,effect_alleles, to harmonise alleles to
  # pval_threshold: only SNPs with P>this threshold will be used to calculate R 
  if(verbose) {
	  cat(
    '    (i) order of column names for effect alleles should be 
    the same as that for estimates and standard errors.
    (ii) it is assumed that dNames is a list of structure
    dNames=list(estNames=c(<x>), seNames=c(<x>), alleleNames=c(<x>))')
  }
  mess=checkBtTerms(merged_data,dNames)
  if(mess!='all good') print(mess)
  Ests=dNames[[1]]
  SEs=dNames[[2]]
  effect_alleles=dNames[[3]]
  merged_data=as.data.frame(merged_data) # not data.table
  e=merged_data[,Ests]
  s=merged_data[,SEs]
  acut=merged_data[,effect_alleles]
  toHarm=acut!=acut[,harmonise_to]
  e[toHarm]=-e[toHarm]
  dcut=e/s
  q=qnorm(1-pval_threshold/2)
  keep=apply(dcut,1,function(h) all(abs(h)<q))
  R=cor(na.omit(dcut[keep,])); colnames(R)=rownames(R)=Ests
  Ncor=sum(keep)
  e=as.matrix(e)
  s=as.matrix(s)
  colnames(s)=colnames(e) # will be required by future functions
  list(R=R,Ncor=Ncor, EstHarm=as.matrix(e), SEHarm=s)
}

#' Plot causal estimates
#'
#' @param fitList list of fitted model objects using one of MRBEE.x()
#' @param exposure_names ordered names of exposures
#' @param qCI quantile for confidence intervals
#' @param binaryOutcome TRUE/FALSE: is outcome binary
#' @param colors colors for causal estimates made using different estimators
#' @param position_to_dodge 0-1 dodging parameter for separation between causal estimates within exposures
#' @keywords 
#' @export
#' @examples
#' causalPlot()

causalPlot=function(fitList,exposure_names,qCI=0.95,binaryOutcome=FALSE,
                    colors=RColorBrewer::brewer.pal(length(fitList),'Spectral'),
                    position_to_dodge=0.8) {
  q=qnorm(qCI+(1-qCI)/2)
  rdf=data.frame()
  meths=c()
  for(i in 1:length(fitList)) {
    method=names(fitList)[i]; meths[i]=method
    Est=c(fitList[[i]]$CausalEstimates)[-1] # take off intercept (always there)
    SE=c(sqrt(diag(fitList[[i]]$VCovCausalEstimates)))[-1] # take off intercept (always)
    rdf=rbind(rdf,data.frame(method,Est,SE,exposure=exposure_names))
  }
  levs=rdf %>% filter(method==meths[1]) %>% arrange(Est) %>% pull(.,exposure)
  rdf$exposure=factor(rdf$exposure,levels=levs)
  if(binaryOutcome) {
    ll=paste0('Odds ratio (', qCI*100, '% CI)')
    p1=ggplot(rdf,aes(x=exposure,y=exp(Est),ymin=exp(Est-q*SE),ymax=exp(Est+q*SE),fill=method)) +
      geom_hline(yintercept=1)
  } else {
    ll=paste0('Causal estimate (', qCI*100, '% CI)')
    p1=ggplot(rdf,aes(x=exposure,y=Est,ymin=Est-q*SE,ymax=Est+q*SE,fill=method)) +
      geom_hline(yintercept=0) +
      labs(y=ll)
  }
  p1=p1+
    geom_pointrange(pch=22,stroke=0.15,position=position_dodge(position_to_dodge),col='gray70') +
    theme_bw() +
    coord_flip() +
    theme(legend.position='bottom') +
    scale_fill_manual('',values=colors) +
    guides(fill=guide_legend(override.aes=list(size=0.8)))
  print(p1)
}

#' Check bias-correction terms can be calculated given data
#'
#' @param merged_data merged_data dataframe of merged exposure(s) and outcome GWAS
#' @param dNames length 3 list of with elements that are vectors of column names in merged_data corresponding to effect sizes, SEs, and effect alleles for all exposures and outcome. Order of elements in list vectors must correspond to the same phenotypes
#' @keywords 
#' @export
#' @examples
#' checkBtTerms()

checkBtTerms=function(merged_data,dNames) {
  merged_data=as.data.frame(merged_data)
  mess='all good'
  cs=sapply(1:3,function(h) class(merged_data[,dNames[[h]]][,1]))
  if(tail(cs,1)!='character') mess=c('effect allele columns named in merged_data should be in final index position of dNames')
  cond1=any(merged_data[,dNames[[1]]]<0)
  cond2=any(merged_data[,dNames[[2]]]<0)
  if(cond1==FALSE | cond2==TRUE) mess=c('effect size columns named in merged_data should be in first index position of dNames and standard errors should be in second index position')
  return(mess)
}

#' Classify outcome loci
#'
#' @param pleioP P-values for horizontal pleiotropy test for outcome
#' @param exposureP P-values for joint test of exposures directly from GWAS
#' @param outcomeP P-values for outcome directly from GWAS
#' @param classifyingLociP P-value to use in classifications of outcome loci into one of direct, mediation, pleiotropy (see published manuscript)
#' @keywords 
#' @export
#' @examples
#' classifyLoci()

classifyLoci=function(pleioP,exposureP,outcomeP,classifyingLociP=5e-8) {
  if(any(c(length(pleioP)==length(exposureP),
           length(pleioP)==length(outcomeP),
           length(exposureP)==length(outcomeP))==FALSE)) stop('not all vectors are of same length')
  if(any(outcomeP>classifyingLociP)) stop('classifyingLociP must be larger than all outcome P-values (outcomeP)')
  m=length(pleioP); type=c()
  for(i in 1:m) {
    if(pleioP[i]<classifyingLociP & exposureP[i]>classifyingLociP & outcomeP[i]<classifyingLociP) type[i]='Direct'
    if(pleioP[i]>classifyingLociP & exposureP[i]<classifyingLociP & outcomeP[i]<classifyingLociP) type[i]='Mediation'
    if(pleioP[i]<classifyingLociP & exposureP[i]<classifyingLociP & outcomeP[i]<classifyingLociP) type[i]='Pleiotropy'
    if(pleioP[i]<classifyingLociP & exposureP[i]<classifyingLociP & outcomeP[i]>classifyingLociP) type[i]='Novel Pleiotropy'
    if(pleioP[i]>classifyingLociP & exposureP[i]>classifyingLociP & outcomeP[i]<classifyingLociP) type[i]='Outcome only'
  }
  return(type)
}

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

#' A function to perform non-iterative MRBEE
#'
#' @param A mxp matrix of standardized GWAS estimates for p exposures and m IVs
#' @param B mxq matrix of standardized GWAS estimates for q outcomes and m IVs
#' @param SigmaUU pxpxm array of measurement error variance-covariance matrices for the p exposures and m IVs.
#' @param SigmaUV pxqxm array of measurement error covariance matrices for the p exposures, q outcomes, and m IVs.
#' @param SigmaVV qxqxm array of measurement error variance-covariance matrices for the q outcomes and m IVs.
#' @param Outliers vector of indices of items to remove from Ahat, Bhat, SigmaUU, SigmaUV, and SigmaVV from their m-index
#' @param perform.pleio logical indicating whether IV-specific horizontal pleiotropy tests should be performed using the statistic Spleio
#' @keywords MRBEE
#' @export
#' @examples
#' EE.estimator()

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

#' Genome-wide horizontal pleiotropy testing
#'
#' @param pD0 output from prepData() ran on full set of exposure and outcome merged GWAS summary statistics
#' @param fit object outputted from one of MRBEE.x()
#' @param verbose TRUE/FALSE should progress be printed
#' @param n_divisions number of times to divide full set of GWAS summary statistics (for computational purposes)
#' @keywords
#' @export
#' @examples
#' genomePleio()

genomePleio=function(pD0,fit,n_divisions=20,verbose=TRUE) {
  # n_divisions: number of parts to divide data into when performing analysis
  # BX=pd$betaX[,,drop=FALSE]
  # BY=pd$betaY[,,drop=FALSE]
  # SX=t(apply(pd$UU,3,function(h) sqrt(c(diag(h)))))
  # SY=sqrt(c(pd$VV)); SY=as.matrix(SY); colnames(SY)=colnames(BY)
  gwRES=matrix(nr=0,nc=8); 
  colnames(gwRES)=c('PleioStat','PleioP','JointExpStat','JointExpP','OutcomeStat','OutcomeP','aHat','aHatHat')
  inds=floor(seq(0,nrow(pD0$betaX),length.out=n_divisions))
  # subR=R[colnames(BX),colnames(BX)]
  if(tail(inds,1)!=nrow(pD0$betaX)) inds[length(inds)]=nrow(pD0$betaX)
  for(i in 1:(length(inds)-1)) {
    di=(inds[i]+1):inds[i+1]
    bx=pD0$betaX[di,,drop=FALSE]
    by=pD0$betaY[di,,drop=FALSE]
    suu=pD0$UU[,,di,drop=FALSE]; # solved=solve(ssum(suu)/nrow(bx)) # much faster
    suv=pD0$UV[,,di,drop=FALSE]
    svv=pD0$VV[,,di,drop=FALSE]
    ostat=c(by)^2/c(svv);ostatp=1-pchisq(ostat,1) # outcome stat and P
    ahathat=cbind(1,bx)%*%fit$CausalEstimates # alphahat and alphathat (alphahat-pleio)
    pdsub=list(bx,by,suu,suv,svv); names(pdsub)=names(pD0)
    sp=Spleio(pdsub,fit)
    # joint test for exposures
    chistats=sapply(1:nrow(bx),function(h) t(bx[h,])%*%solve(suu[,,h])%*%bx[h,])
    gwRES=rbind(gwRES,cbind(sp$Stats,sp$Ps, chistats,1-pchisq(chistats,ncol(bx)),ostat,ostatp,c(by),c(ahathat)))
    if(verbose) {
      if(i%in%(1:length(inds))) print(paste0(round(i/length(inds)*100),'% complete'))
    }
  }
  as.data.frame(gwRES)
}

#' Plotting results from genome-wide horizontal pleiotropy testing
#'
#' @param gwSp object output of genomePleio()
#' @param leadSNPs rsIDs of lead outcome SNPs
#' @param CHR_BP vector of type "CHR_BP" positions corresponding to SNPs used in genome-wide horizontal pleiotropy testing
#' @param GWsubInds indices in full merged GWAS data to use in plotting (only change if you want plots for only a subset of all available SNPs)
#' @param classifyingLociP P-value to use in classifications of outcome loci into one of direct, mediation, pleiotropy (see published manuscript)
#' @param manhattanPUpper only SNPs with P<this threshold will be plotted in the Manhattan plot
#' @keywords 
#' @export
#' @examples
#' gwPleioPlot()

gwPleioPlot=function(gwSp,leadSNPs,CHR_BP,GWsubInds=1:length(gwSp$PleioStat),
                     classifyingLociP=5e-8,manhattanPUpper=0.05) {
  if(length(CHR_BP)!=length(gwSp$PleioStat)) stop('length of CHR_BP vector is not of same length as number of SNPs tested in gwSp')
  if(is.character(leadSNPs)) stop('leadSNPs should be a vector of row indices (integers) corresponding to SNPs tested in gwSp')
  ### QQ plot
  # chisq joint test for exposures is present in gwSp
  p=0; error=1; # finding degrees of freedom it must have used
  while(error>0.000001) { p=p+1; error=gwSp$JointExpStat[1]-qchisq(1-gwSp$JointExpP[1],p) }
  Stats=data.frame(pleio=c(gwSp$PleioStat),exposures=c(gwSp$JointExpStat),outcome=c(gwSp$OutcomeStat))
  Stats=Stats[GWsubInds,]
  Stats=na.omit(Stats)
  Stats$expectedExposures=qchisq(ppoints(nrow(Stats)),p) # already sorted
  Stats$expectedOutcomePleio=qchisq(ppoints(nrow(Stats)),1) # already sorted
  qq1=ggplot(Stats,aes(expectedExposures,sort(exposures),color='a')) +
    geom_abline(intercept=0,slope=1) +
    geom_point(size=0.75,) +
    labs(x=bquote('expected '*chi^2*'('*.(p)*')'),
         y=bquote('observed '*chi^2*'('*.(p)*')')) +
    theme_bw() +
    scale_color_manual('',labels=c('Joint test for exposures'),values='indianred') +
    theme(legend.text=element_text(size=12),legend.position='bottom') +
    guides(color=guide_legend(override.aes=list(size=3)))
  qq2=ggplot(Stats,aes(expectedOutcomePleio,sort(outcome),color='a')) +
    geom_abline(intercept=0,slope=1) +
    geom_point(size=0.75) +
    geom_point(aes(expectedOutcomePleio,sort(pleio),color='b'),size=0.75,data=Stats) +
    labs(x=bquote('expected '*chi^2*'(1)'),
         y=bquote('observed '*chi^2*'(1)')) +
    theme_bw() +
    scale_color_manual('',labels=c('Outcome GWAS','Horizontal pleiotropy'),
                       values=c('black','royalblue')) +
    theme(legend.text=element_text(size=12),legend.position='bottom') +
    guides(color=guide_legend(override.aes=list(size=3)))
  arrangedqq=ggpubr::ggarrange(qq2,qq1,nrow=1,ncol=2)
  
  ### locus plot
  df=data.frame(obsPleio=gwSp$PleioP,obsY=gwSp$OutcomeP,obsX=gwSp$JointExpP,
                ahat=gwSp$aHat,ahatHat=gwSp$aHatHat)
  dfsub=df[leadSNPs,]
  dfsub$type=classifyLoci(dfsub$obsPleio,dfsub$obsX,dfsub$obsY,classifyingLociP=classifyingLociP)
  rho=cor.test(dfsub$ahat[dfsub$type=="Mediation"],dfsub$ahatHat[dfsub$type=="Mediation"],method='spearman')
  rhoP=rho$p.value; rho=rho$estimate
  lab1=expression(underline('Spearman '*rho))
  rho2=cor.test(dfsub$ahat[dfsub$type!="Mediation"],dfsub$ahatHat[dfsub$type!="Mediation"])
  rhoP2=rho2$p.value; rho2=rho2$estimate
  rhoP=ifelse(rhoP<0.001,'<0.001',paste0('=',as.character(round(rhoP,3))))
  rhoP2=ifelse(rhoP2<0.001,'<0.001',paste0('=',as.character(round(rhoP2,3))))
  rho1=as.character(round(rho,2))
  rho2=as.character(round(rho2,2))
  lab2=paste0('Mediation: ', rho1, ', P',rhoP)
  lab3=paste0('Non-mediation: ', rho2, ', P',rhoP2)
  fit=lm(ahat~ahatHat,data=dfsub%>%filter(type=="Mediation"))
  p3=ggplot(dfsub,aes(ahatHat,ahat,fill=type)) +
    geom_abline(intercept=0,slope=1,linetype="dotted") +
    geom_hline(yintercept=0,linetype="dotted") +
    geom_vline(xintercept=0,linetype="dotted") +
    geom_point(pch=21,size=2) +
    theme_bw() +
    # scale_fill_brewer(palette="Accent") +
    scale_fill_brewer(palette="Set2") +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(override.aes=list(size=3),title="")) +
    theme(legend.text=element_text(size=12)) +
    labs(x='total SNP effect through exposures',y='total SNP effect on outcome')
  bb=ggplot_build(p3)
  medcol=unique(bb$data[[4]]["fill"]$fill[dfsub$type=="Mediation"]) # color of mediation points
  labloc=c(min(dfsub$ahatHat),max(dfsub$ahat))*c(0.75,0.95)
  fulllab=c(lab1, '\n', lab2, '\n', lab3)
  defaultW=getOption("warn");options(warn=-1)  
  p3=p3+geom_abline(intercept=coef(fit)[1],slope=coef(fit)[2],color=medcol) +
    annotate('text',x=labloc[1],y=labloc[2],label=lab1,size=3.5) +
    annotate('text',x=labloc[1],y=labloc[2]*0.9,label=lab2,size=3.5) +
    annotate('text',x=labloc[1],y=labloc[2]*0.8,label=lab3,size=3.5)
  options(warn=defaultW)
  
  ### manhattan plots (top and bottom like in MRBEE paper)
  res=strsplit(CHR_BP,"_",fixed=TRUE)
  m=t(matrix(t(unlist(res)),nr=2,nc=length(CHR_BP)))
  m=apply(m,2,as.numeric);colnames(m)=c('CHR','BP') # faster than sapply()
  dfman=df %>% mutate(CHR=m[,1],BP=m[,2])
  dfman=dfman %>% arrange(CHR,BP)
  dfman=dfman[GWsubInds,]
  dfman=dfman %>% filter(obsY<manhattanPUpper | obsPleio<manhattanPUpper)
  # for code below, see: https://r-graph-gallery.com/101_Manhattan_plot.html
  don=dfman %>% group_by(CHR) %>% summarise(chr_len=max(BP,na.rm=T)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(.,dfman,by="CHR") %>%
    arrange(CHR,BP) %>%
    mutate(BPcum=BP+tot)
  axisdf=don %>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)
  # top: outcome GWAS
  ptop=ggplot(don,aes(BPcum,-log10(obsY))) +
    geom_point( aes(color=as.factor(CHR)), alpha=1, size=1.3) +
    scale_color_manual(values=rep(c("skyblue", "royalblue"), 22)) +
    scale_x_continuous(label=axisdf$CHR,breaks=axisdf$center) +
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x=NULL,y=expression('-log'[10]*'(P-value)'),title='Outcome GWAS') +
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size=7),plot.title=element_text(hjust=0.5)) +
    geom_hline(yintercept=-log10(5e-5),color='red') +
    geom_hline(yintercept=-log10(5e-8),color='black')
  # bottom: pleiotropy testing
  pbottom=ggplot(don,aes(BPcum,-log10(obsPleio))) +
    geom_point(aes(color=as.factor(CHR)), alpha=1, size=1.3) +
    scale_color_manual(values=rep(c("skyblue", "royalblue"), 22)) +
    scale_x_continuous(label=axisdf$CHR,breaks=axisdf$center) +
    theme_bw() +
    scale_y_reverse() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(x=NULL,y=expression('-log'[10]*'(P-value)'),title=NULL,caption='Horizontal pleiotropy test') +
    theme(plot.caption = element_text(hjust=0.5, size=rel(1.2))) +
    theme(axis.text.x=element_blank()) +
    geom_hline(yintercept=-log10(5e-5),color='red') +
    geom_hline(yintercept=-log10(5e-8),color='black')
  arranged=ggpubr::ggarrange(ptop,pbottom,nrow=2)
  ### done
  outList=list(QQ=arrangedqq,Manhattan=arranged,LocusPlot=p3)
  return(outList)
}

#' Helper
#'
#' @param .array 
#' @keywords
#' @export
#' @examples
#' array.to.list()
array.to.list=function(.array) lapply(seq(dim(.array)[3]), function(x) .array[ , , x])

#' Helper
#'
#' @param R
#' @param x
#' @keywords
#' @export
#' @examples
#' ff()
ff=function(R,x) sapply(x,function(h) which(h==colnames(R)))

#' Helper
#'
#' @param C
#' @param Est
#' @param Var
#' @keywords
#' @export
#' @examples
#' fQ()
fQ=function(C, Est, Var) t(C %*% c(Est)) %*% solve(C %*% Var %*% t(C)) %*% (C %*% c(Est))

#' Helper
#'
#' @param Bhati
#' @param Thetahat
#' @param Ahati
#' @param Upsilonhati
#' @keywords
#' @export
#' @examples
#' fS()
fS=function(Bhati, Thetahat, Ahati, Upsilonhati) t(Bhati-t(Thetahat) %*% Ahati) %*% solve(Upsilonhati) %*% (Bhati-t(Thetahat) %*% Ahati)

#' Helper
#'
#' @param .list
#' @keywords
#' @export
#' @examples
#' ssum()
ssum=function(.list) {
  if(is.array(.list)) .list=array.to.list(.list)
  return(Reduce("+", .list))
}

#' Helper
#'
#' @param A
#' @param B
#' @param SigmaUU
#' @param SigmaUV
#' @param SigmaVV
#' @param Outliers
#' @keywords
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

#' Helper
#'
#' @param A
#' @keywords
#' @export
#' @examples
#' vec()
vec=function(A){
  a=as.vector(as.matrix(A))
  return(a)
}

#' Helper
#'
#' @param A
#' @keywords
#' @export
#' @examples
#' demean()
demean=function(A){
  a=colMeans(A)
  b=t(t(A)-a)
  return(b)
}

#' Helper
#'
#' @param A
#' @param B
#' @keywords
#' @export
#' @examples
#' cdiag()
cdiag=function(A,B){
  # p1=ncol(asA);p2=ncol(B)
  p1=ncol(as.matrix(A));p2=ncol(as.matrix(B))
  if(sum(p1)!=0&sum(p2)!=0){
    C=diag(p1+p2)
    C[1:p1,1:p1]=A
    C[(p1+1):(p1+p2),(p1+1):(p1+p2)]=B
  }
  if(sum(p1)!=0&sum(p2)==0){
    C=diag(p1+1)
    C[1:p1,1:p1]=A
    C[p1+1,p1+1]=B
  }
  if(sum(p1)==0&sum(p2)!=0){
    #C=diag(p1+1)
    #C[1,1]=A
    #C[2:(p2+1),2:(p2+1)]=B  
    C=matrix(0,p1+p2,p1+p2)
    C[2:nrow(C),2:ncol(C)]=B
  }
  if(sum(p1)==0&sum(p2)==0){
    C=diag(2);C[1,1]=A;C[2,2]=B
  }
  return(C)
}

#' Helper
#'
#' @param Vbb
#' @param Vuu
#' @param Ruv
#' @param theta
#' @param hy
#' @keywords
#' @export
#' @examples
#' mvvarv()
mvvarv=function(Vbb,Vuu,Ruv,theta,hy){
  Vu=diag(Vuu)
  b=2*sum(theta*Ruv*sqrt(Vu))
  c=t(theta)%*%(Vbb+Vuu)%*%theta-t(theta)%*%Vbb%*%theta/hy
  vv=(-b+sqrt(b^2-4*c))/2
  return(vv^2)
}

#' Helper
#'
#' @param n
#' @param mean
#' @param sigma
#' @keywords
#' @export
#' @examples
#' sigma2normal()
sigma2normal=function(n,mean,sigma){
  r=rnorm(5*n,mean,sigma)
  ind=which(abs(r)<sigma)
  r=r[-ind]
  if(length(r)<n){
    r1=rnorm(5*n,mean,sigma)
    ind1=which(abs(r1)<sigma)
    r1=r1[-ind]
    r=c(r,r1)
  }
  m=length(r)
  r=r[sample(c(1:m),n,replace=F)]
  return(r)
}

#' Helper
#'
#' @param vbb
#' @param vuu
#' @param ruv
#' @param theta
#' @param hy
#' @keywords
#' @export
#' @examples
#' varv()
varv=function(vbb,vuu,ruv,theta,hy){
  b=2*theta*ruv*sqrt(vuu)
  c=(vbb+vuu)*theta^2-vbb*theta^2/hy
  vv=(-b+sqrt(b^2-4*c))/2
  return(vv^2)
}

#' Helper
#'
#' @param H
#' @keywords
#' @export
#' @examples
#' ginv()
ginv=function(H){
  h=eigen(H)
  d=h$values
  v=h$vectors
  ind=which(d>0)
  d[ind]=1/d[ind]
  d[-ind]=0
  H=v%*%diag(d)%*%t(v)
  return(H)
}

#' Helper
#'
#' @param b
#' @keywords
#' @export
#' @examples
#' svar()
svar=function(b){
  b=b-mean(b)
  y=sum(b^2)
  return(y)
}

#' Helper
#'
#' @param a
#' @keywords
#' @export
#' @examples
#' diags()
diags=function(a){
  if(length(a)>1) A=diag(a)
  if(length(a)==1) A=a
  return(A)
}

#' Helper
#'
#' @param a
#' @param ase
#' @param a0
#' @keywords
#' @export
#' @examples
#' iscov()
iscov=function(a,ase,a0){
  d1=abs(a-a0)
  b=as.numeric(d1<(2*ase))
  return(b)
}

#' Helper
#'
#' @param a
#' @keywords
#' @export
#' @examples
#' robustse()
robustse=function(a){
  b=1.483*median(abs(a-median(a)))
  return(b)
}

#' Helper
#'
#' @param a
#' @param lam
#' @param gamma
#' @keywords
#' @export
#' @examples
#' dSCAD()
dSCAD=function(a,lam,gamma=3.7){
  a=abs(a)
  z=a
  z[a<lam]=lam
  z[a>lam]=(gamma*lam-z[a>lam])/(gamma-1)
  z[a>(gamma*lam)]=0  
  return(z)
}

#' Helper
#'
#' @param a
#' @param lam
#' @param gamma
#' @keywords
#' @export
#' @examples
#' dMCP()
dMCP=function(a,lam,gamma=3){
  a=abs(a)
  z=lam-a/gamma
  z[a>(gamma*lam)]=0  
  return(z)
}

#' Helper
#'
#' @param a
#' @param lam
#' @param gamma
#' @keywords
#' @export
#' @examples
#' tSCAD()
tSCAD=function(a,lam,gamma=3.7){
  b=abs(a)
  c=b
  ind=which(abs(b)<=(gamma*lam)&abs(b)>(2*lam))
  c[ind]=soft(b[ind],gamma*lam/(gamma-1))/(1-1/(gamma-1))
  ind=which(abs(b)<=(2*lam))
  c[ind]=soft(b[ind],lam)
  return(c*sign(a))
}

#' Helper
#'
#' @param a
#' @param lam
#' @param gamma
#' @keywords
#' @export
#' @examples
#' tMCP()
tMCP=function(a,lam,gamma=3){
  b=abs(a)
  c=soft(b,lam)/(1-1/gamma)
  c[b>(gamma*lam)]=b[b>(gamma*lam)]
  c=c*sign(a)
  return(c)
}

#' Helper
#'
#' @param a
#' @param lam
#' @keywords
#' @export
#' @examples
#' tlasso()
tlasso=function(a,lam){
  b=soft(a,rep(lam,length(a)))
  return(b)
}

#' Helper
#'
#' @param a
#' @param b
#' @keywords
#' @export
#' @examples
#' soft()
soft=function(a,b){
  c=abs(a)-b
  c[c<0]=0
  c=c*sign(a)
  return(c)
}

#' Helper
#'
#' @param x
#' @param theta
#' @param Rww
#' @param rwy
#' @param q
#' @param method
#' @keywords
#' @export
#' @examples
#' mlqeweight()
mlqeweight=function(x,theta,Rww,rwy,q=0.95,method="ordinal"){
  p=length(Rww[,1])
  Rwxy=cdiag(Rww,1)
  Rwxy[p+1,1:p]=Rwxy[1:p,p+1]=rwy # changed rxy to rwy
  theta=c(theta,-1)
  varx=t(theta)%*%Rwxy%*%theta
  if(method=="robust") varx=(1.483*median(abs(x-median(x))))^2
  f=dnorm(x,mean=0,sd=sqrt(varx))
  weight=as.vector(f^(1-q))
  weight=as.vector(weight/sum(weight))
  return(weight)
}

#' Helper
#'
#' @param x
#' @param theta
#' @param Rww
#' @param rwy
#' @param se.est
#' @param FDR
#' @param adjust.method
#' @keywords
#' @export
#' @examples
#' imrpdetect()
imrpdetect=function(x,theta,Rww,rwy,se.est="ordinal",FDR=T,adjust.method="Sidak"){
  p=length(Rww[1,])
  if(se.est=="robust"){ 
    varx=(1.483*median(abs(x-median(x))))^2
  }else{
    Rwxy=cdiag(Rww,1)
    Rwxy[p+1,1:p]=Rwxy[1:p,p+1]=rxy
    theta=c(theta,-1)
    varx=as.numeric(t(theta)%*%Rwxy%*%theta)
  }
  pv=1-pchisq(x^2/varx,1)
  if(FDR==T){
    pv=FDRestimation::p.fdr(pvalues=pv,adjust.method=adjust.method)$fdrs
  }
  return(as.vector(pv))
}

#' Helper
#'
#' @param by
#' @param bW
#' @param theta1
#' @param theta1
#' @param ind1
#' @param ind2
#' @param outlier
#' @param thres
#' @keywords
#' @export
#' @examples
#' normvoting()
normvoting=function(by,bW,theta1,theta2,ind1,ind2,outlier=F,thres=0.05){
  y1=by[ind1];y2=by[ind2];W1=bW[ind1,];W2=bW[ind2,]
  e1=y1-W1%*%theta1;e2=y2-W2%*%theta2;
  v1=robustse(e1);v2=robustse(e2);
  r1=by-bW%*%theta1;r2=by-bW%*%theta2;
  p1=dnorm(r1,sd=v1);p2=dnorm(r2,sd=v2);
  p3=p1+p2
  p1=p1/p3;p2=p2/p3
  p3=p1-p2
  ind1=which(p3>=0)
  ind2=which(p3<0)
  
  if(length(ind1)<length(ind2)){
    ind3=ind1;ind1=ind2;ind2=ind3
    p3=p1;p1=p2;p2=p3
  }
  
  p1=sqrt(p1/sum(p1))
  p2=sqrt(p2/sum(p2))
  A=list()
  A$ind1=c(ind1)
  A$ind2=c(ind2)
  A$p1=c(p1)
  A$p2=c(p2)
  
  if(outlier==T){
    pv1=1-pchisq(r1^2/v1^2,1);pv2=1-pchisq(r2^2/v2^2,1)
    pv=vecmax(pv1,pv2)
    pv=FDRestimation::p.fdr(pv)$fdrs
    indo=pv>thres
    A$p1=A$p1*indo
    A$p2=A$p2*indo
    A$outlier=which(indo==0)
  }
  return(A)
}

#' Helper
#'
#' @param a
#' @param b
#' @keywords
#' @export
#' @examples
#' vecmax()
vecmax=function(a,b){
  n=length(a)
  d=a
  for(i in 1:n){
    d[i]=max(a[i],b[i])
  }
  return(d)
}

#' Helper
#'
#' @param A
#' @keywords
#' @export
#' @examples
#' coldSD()
colSD=function(A){
  B=A
  p=length(A[1,])
  for(i in 1:p){
    B[,i]=robustse(A[,i])  
  }
  return(B)
}

#' Helper
#'
#' @param by
#' @param bX
#' @param theta
#' @param Suu
#' @param Suv
#' @param h
#' @param tau
#' @keywords
#' @export
#' @examples
#' smooth.quantile()
smooth.quantile=function(by,bX,theta,Suu,Suv,h,tau=0.5){
  e=as.vector(by-bX%*%theta)
  sigma=c(t(theta)%*%Suu%*%theta-2*sum(Suv*theta))
  resj=sapply(1:length(by), function(i) {
    int_=integrate(int.median,lower=0,upper=1/h,e=e[i],sigma=sigma,h=h)$value
    (by[i]-bX[i,]%*%theta)*(tau-1/2)+1/pi*int_
  })
  sum(resj)
}

#' Helper
#'
#' @param theta
#' @param by
#' @param bX
#' @param Suu
#' @param Suv
#' @param h
#' @param tau
#' @keywords
#' @export
#' @examples
#' gr.smooth.quantile()
gr.smooth.quantile=function(theta,by,bX,Suu,Suv,h,tau=0.5){
  e=as.vector(by-bX%*%theta)
  sigma=c(t(theta)%*%Suu%*%theta-2*sum(Suv*theta))
  W=quantileweight(e,sigma,h)
  
  G=-c((tau-0.5)*colSums(bX))-c(t(bX)%*%W[,1])+sum(W[,2])*2*c(Suu%*%theta-Suv)
  return(G)
}

#' Helper
#'
#' @param y
#' @param e
#' @param sigma
#' @keywords
#' @export
#' @examples
#' weightfun1()
weightfun1=function(y,e,sigma){
  a=sin(y*e)/y+e*cos(y*e)+y*sigma*sin(y*e)  
  b=exp(y^2*sigma^2/2)  
  return(a*b)  
}

#' Helper
#'
#' @param y
#' @param e
#' @param sigma
#' @keywords
#' @export
#' @examples
#' weightfun2()
weightfun2=function(y,e,sigma){
  a=-cos(y*e)*exp(y^2*sigma^2/2)/2
  b=y^2/2*(1/y*e*sin(y*e)-sigma*cos(y*e))*exp(y^2*sigma/2)
  return(a+b)
}

#' Helper
#'
#' @param e
#' @param sigma
#' @param h
#' @keywords
#' @export
#' @examples
#' quantileweight()
quantileweight=function(e,sigma,h){
  w1=sapply(1:length(e),function(i){
    integrate(weightfun1,lower=0,upper=1/h,e=e[i],sigma=sigma)$value
  })
  w2=sapply(1:length(e),function(i){
    integrate(weightfun2,lower=0,upper=1/h,e=e[i],sigma=sigma)$value
  })
  W=cbind(w1,w2)
  return(W)
}

#' Helper
#'
#' @param by
#' @param bX
#' @param theta
#' @param tau
#' @param Suu
#' @param Suv
#' @param h
#' @keywords
#' @export
#' @examples
#' medianerror()
medianerror=function(by,bX,theta,tau,Suu,Suv,h){
  e=c(by-bX%*%theta)
  sigma=t(theta)%*%Suu%*%theta-2*t(theta)%*%Suv
  W=quantileweight(e,sigma,h)
  w1=c(W[,1])
  w2=c(W[,2])
  E=-bX*(tau-0.5)-bX*w1+kronecker(matrix(1,length(e),1),2*t(as.vector(Suv-Suu%*%theta)))*w2
  return(E)
}

#' Helper
#'
#' @param g
#' @param e
#' @param sigma
#' @param h
#' @keywords
#' @export
#' @examples
#' int.median()
int.median=function(g,e,sigma,h){
  s=c(1/g*e*sin(g*e)-sigma*cos(g*e))*exp(g^2*sigma/2)
  return(s)
}

#' Internal MRBEE-IMRP function
#'
#' @param by
#' @param bX
#' @param Rxx
#' @param rxy
#' @param max.iter
#' @param max.eps
#' @param pv.thres
#' @param se.est
#' @param FDR
#' @param adjust.method
#' @param intercept
#' @keywords
#' @export
#' @examples
#' imrp.mrbee.internal()

imrp.mrbee.internal=function(by,bX,Rxx,rxy,max.iter=15,max.eps=1e-3,pv.thres=0.05,
                             se.est="robust",FDR=T,adjust.method="Sidak",intercept=T){
  n=length(by)
  bW=cbind(1,bX)
  Rww=cdiag(0,Rxx)
  rwy=c(0,rxy)
  Hinv=solve(t(bW)%*%bW-n*Rww)
  theta=Hinv%*%(t(bW)%*%by-n*rwy)
  theta1=theta*0
  error=sqrt(sum((theta-theta1)^2))
  iter=0
  while(error>max.eps&iter<max.iter){
    theta1=theta
    e=vec(by-bW%*%theta)  
    pv=imrpdetect(x=e,theta=theta,Rww=Rww,rwy=rwy,se.est=se.est,FDR=FDR,adjust.method=adjust.method)
    ind=which(pv>pv.thres)
    Hinv=solve(t(bW[ind,])%*%bW[ind,]-length(ind)*Rww)
    theta=Hinv%*%(t(bW[ind,])%*%by[ind]-length(ind)*rwy)
    theta[1]=theta[1]*intercept
    iter=iter+1
    error=sqrt(sum((theta-theta1)^2))
  }
  E=-bW[ind,]*e[ind]+kronecker(matrix(1,length(ind),1),t(-vec(Rww%*%theta)+rwy))
  V=t(E)%*%E
  covtheta=Hinv%*%V%*%Hinv
  A=list()
  A$theta=theta
  A$covtheta=covtheta
  
  r=by-bW%*%theta
  r[ind]=0
  A$delta=r
  return(A)
}

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

#' MRBEE with IMRP adjustment
#'
#' @param prepDataList output of prepData() on the IV set
#' @param PleioPThreshold P-value threshold for determining if a specific IV has sufficient evidence of horizontal pleiotropy to be removed from causal estimation
#' @param FDR TRUE/FALSE should Benjamini-Hochberg FDR adjustment be used when identifying horizontally pleiotropic IVs
#' @param FDR.alpha Type I error rate for FDR adjustment
#' @param inds indices ot IV set to use - only change if you want to use an a-priori defined subset of all IVs
#' @param UU pxpxm array of measurement error variance-covariance matrices for the p exposures and m IVs.
#' @param UV pxqxm array of measurement error covariance matrices for the p exposures, q outcomes, and m IVs.
#' @param VV qxqxm array of measurement error variance-covariance matrices for the q outcomes and m IVs.
#' @keywords 
#' @export 
#' @examples
#' MRBEE.IMRP()

MRBEE.IMRP=function(prepDataList, PleioPThreshold=0.05, FDR=FALSE, FDR.alpha=0.05, 
                    inds=1:nrow(prepDataList$betaX)) {
  bx=prepDataList$betaX[inds,,drop=FALSE];
  by=prepDataList$betaY[inds,,drop=FALSE];
  UU=prepDataList$UU[,,inds,drop=FALSE];
  UV=prepDataList$UV[,,inds,drop=FALSE];
  VV=prepDataList$VV[,,inds,drop=FALSE]
  it=iter.estimator(A=bx,B=by,SigmaUU=UU,SigmaUV=UV,SigmaVV=VV,PleioPThreshold=PleioPThreshold,FDR=FDR,FDR.alpha=FDR.alpha,Outliers=NULL)
  it$method="imrp"
  return(it)
}

#' MRBEE with IPOD adjustment
#'
#' @param pd output of prepData() on the IV set
#' @param max.iter maximum number of iterations to perform until convergence
#' @param max.eps tolerance (smaller values are stricter)
#' @param penalty type of penalty, one of "SCAD", "MCP", "lasso"
#' @param tau tuning parameter lambda for SCAD, MCP, and lasso penalties
#' @param intercept TRUE/FALSE should an intercept be included in causal estimation    
#' @keywords 
#' @export 
#' @examples
#' MRBEE.IPOD()

MRBEE.IPOD=function(pd,max.iter=15,max.eps=1e-3,penalty="SCAD",tau=3,intercept=T){
  by=pd$betaY
  bX=pd$betaX
  Rxx=ssum(pd$UU)/length(pd$betaY)
  rxy=ssum(pd$UV)/length(pd$betaY)
  n=length(by)
  bW=cbind(1,bX)
  Rww=cdiag(0,Rxx)
  rwy=c(0,rxy)
  Hinv=solve(t(bW)%*%bW-n*Rww)
  theta=Hinv%*%(t(bW)%*%by-n*rwy)
  n=length(by)
  p=length(bW[1,])
  r=vec(by-bW%*%theta)
  ################### SCAD ####################
  if(penalty=="SCAD"){
    theta1=theta*0
    error=sqrt(sum((theta1-theta)^2))/sqrt(p)
    iter=0
    while(error>max.eps&iter<max.iter){
      theta1=theta
      r=vec(by-bW%*%theta)
      delta=tSCAD(r,tau) 
      ind=(delta==0) 
      theta=vec(solve(t(bW)%*%bW-sum(ind)*Rww)%*%(t(bW)%*%(by-delta)-sum(ind)*rwy))
      theta[1]=theta[1]*intercept
      error=sqrt(sum((theta1-theta)^2))/sqrt(p)
      iter=iter+1
    }
  }  
  
  if(penalty=="MCP"){
    theta1=theta*0
    error=sqrt(sum((theta1-theta)^2))/sqrt(p)
    iter=0
    while(error>max.eps&iter<max.iter){
      theta1=theta
      r=vec(by-bW%*%theta)
      delta=tMCP(r,tau) 
      ind=(delta==0) 
      theta=vec(solve(t(bW)%*%bW-sum(ind)*Rww)%*%(t(bW)%*%(by-delta)-sum(ind)*rwy))
      theta[1]=theta[1]*intercept
      error=sqrt(sum((theta1-theta)^2))/sqrt(p)
      iter=iter+1
    }
  }  
  
  if(penalty=="lasso"){
    theta1=theta*0
    error=sqrt(sum((theta1-theta)^2))/sqrt(p)
    iter=0
    while(error>max.eps&iter<max.iter){
      theta1=theta
      r=vec(by-bW%*%theta)
      delta=tlasso(r,tau) 
      ind=(delta==0) 
      theta=vec(solve(t(bW)%*%bW-sum(ind)*Rww)%*%(t(bW)%*%(by-delta)-sum(ind)*rwy))
      theta[1]=theta[1]*intercept
      error=sqrt(sum((theta1-theta)^2))/sqrt(p)
      iter=iter+1
    }
  }  
  ################### inference #################
  D=diag(n)
  ind=which(delta!=0)
  r=vec(by-bW%*%theta-delta)
  if(length(ind)>0){
    D=D[,ind]
    bZ=cbind(bW,D)
    thetadelta=c(theta,delta[which(delta!=0)])
    Rzz=cdiag(Rww,diag(ind)*0)
    rzy=c(rwy,0*ind)
    E=-bZ*r
    E[-ind,]=E[-ind,]+kronecker(matrix(1,n-length(ind),1),t(vec(rzy-Rzz%*%thetadelta)))
    V=t(E)%*%E
    H=t(bZ)%*%bZ-(n-length(ind))*Rzz
    S=solve(H)%*%V%*%solve(H)
    S=S[1:p,1:p]
  }else{
    E=-bW*r+kronecker(matrix(1,n,1),t(vec(rwy-Rww%*%theta))) 
    H=t(bW)%*%bW-n*Rww
    V=t(E)%*%E
    S=solve(H)%*%V%*%solve(H)
  }
  
  A=list()
  A$CausalEstimates=theta
  A$delta=delta
  A$VCovCausalEstimates=S
  A$method='ipod'
  return(A)
}

#' MRBEE with quantile estimation
#'
#' @param pd output of prepData() on the IV set
#' @param tau quantile to estimate (tau=0.5 is median regression)
#' @param h tuning parameter in (0,Inf), should be kept relatively large (> ~0.5). Computation is faster for larger values
#' @param intercept TRUE/FALSE should an intercept be included in causal estimation    
#' @param gradient TRUE/FALSE should the Hessian of the loss function be used in optimization
#' @keywords 
#' @export 
#' @examples
#' MRBEE.Median()

MRBEE.Median=function(pd,tau=0.5,h=1,intercept=T,gradient=F){
  by=pd$betaY
  bX=pd$betaX
  Rxx=ssum(pd$UU)/length(pd$betaY)
  rxy=ssum(pd$UV)/length(pd$betaY)
  n=length(by)
  bW=cbind(1,bX)
  Rww=cdiag(0,Rxx)
  rwy=c(0,rxy)
  if(intercept==T){
    theta=MRBEE.IMRP(pd,PleioPThreshold=0.05/nrow(pd$betaX),FDR=FALSE)$CausalEstimates
    # theta=imrp.mrbee(by,bX,Rxx,rxy)$theta
    # medianrun=optim(par=theta,smooth.quantile,by=by,bX=bW,Suu=Rww,Suv=rwy,h=h,tau=tau,hessian=TRUE)
    if(gradient) {
      medianrun=optim(par=theta,fn=smooth.quantile,by=by,bX=bW,Suu=Rww,Suv=rwy,h=h,tau=tau,hessian=TRUE)
    } else {
      medianrun=optim(par=theta,fn=smooth.quantile,gr=gr.smooth.quantile,by=by,bX=bW,Suu=Rww,Suv=rwy,h=h,tau=tau,hessian=TRUE)
    }
    theta=medianrun$par
    H=medianrun$hessian
    #E=medianerror(by=by,bX=bW,theta=theta,tau=tau,Suu=Rww,Suv=rwy,h=h)
  }else{
    theta=imrp.mrbee(by,bX,Rxx,rxy)$theta[-1] 
    # medianrun=optim(par=theta,smooth.quantile,by=by,bX=bX,Suu=Rxx,Suv=rxy,h=h,tau=tau,hessian=TRUE)
    if(gradient) {
      medianrun=optim(par=theta,fn=smooth.quantile,gr=gr.smooth.quantile,by=by,bX=bX,Suu=Rxx,Suv=rxy,h=h,tau=tau,hessian=TRUE)
    } else{
      medianrun=optim(par=theta,fn=smooth.quantile,by=by,bX=bX,Suu=Rxx,Suv=rxy,h=h,tau=tau,hessian=TRUE)
    }
    theta=medianrun$par
    #H=medianrun$hessian
    #E=medianerror(by=by,bX=bx,theta=theta,tau=tau,Suu=Rxx,Suv=rxy,h=h)
  }
  #V=t(E)%*%E
  thetacov=solve(H)
  #%*%V%*%solve(H)
  A=list()
  A$CausalEstimates=theta
  A$VCovCausalEstimates=thetacov
  A$method='median'
  return(A)
}

#' MRBEE assuming a mixture distribution for IVs
#'
#' @param pd output of prepData() on the IV set
#' @param max.iter maximum number of iterations to perform until convergence
#' @param max.eps tolerance (smaller values are stricter)
#' @param FDR TRUE/FALSE should Benjamini-Hochberg FDR adjustment be used when identifying horizontally pleiotropic IVs
#' @param intercept1 TRUE/FALSE should an intercept be included in causal estimation for the first of 2 IV sets
#' @param intercept2 TRUE/FALSE should an intercept be included in causal estimation for the second of 2 IV sets
#' @param bic TRUE/FALSE should BIC be calculated
#' @param aic TRUE/FALSE should AIC be calculated  
#' @keywords 
#' @export 
#' @examples
#' MRBEE.Mix()

MRBEE.Mix=function(pd,max.iter=15,max.eps=1e-3,FDR=T,intercept1=T,intercept2=F,bic=T,aic=F){
  by=pd$betaY
  bX=pd$betaX
  Rxx=ssum(pd$UU)/length(pd$betaY)
  rxy=ssum(pd$UV)/length(pd$betaY)
  n=length(by)
  bX=as.matrix(bX)
  bW=cbind(1,bX)
  p=length(bW[1,])
  Rww=cdiag(0,Rxx)
  rwy=c(0,rxy)
  
  fit1=imrp.mrbee.internal(by,bX,Rxx,rxy,FDR=F,intercept=intercept1)
  delta1=fit1$delta
  theta1=fit1$theta
  ind1=which(delta1==0)
  ind2=which(delta1!=0)
  
  # fit2=imrp.mrbee.internal(pd,PleioPThreshold=0,inds=ind2)
  fit2=imrp.mrbee.internal(by[ind2],bX[ind2,],Rxx,rxy,intercept=intercept2)
  theta2=fit2$theta
  vote=normvoting(by=by,bW=bW,theta1=theta1,theta2=theta2,ind1=ind1,ind2=ind2)
  ind1=vote$ind1
  ind2=vote$ind2
  p1=vote$p1
  p2=vote$p2
  theta11=theta1*0
  error=sqrt(sum(theta11-theta1)^2)
  iter=0
  
  while(error>max.eps & iter<max.iter){
    theta11=theta1
    # fit1=weight.mrbee.internal(by,bX,weight=p1,Rxx,rxy,intercept=intercept1)
    fit1=weight.mrbee.internal(pd,weight=p1,intercept=intercept1)
    # theta1=fit1$theta
    theta1=fit1$CausalEstimates
    if(sum(ind2)>0){
      # fit2=weight.mrbee.internal(by,bX,weight=p2,Rxx,rxy,intercept=intercept2)
      fit2=weight.mrbee.internal(pd,weight=p2,intercept=intercept2)
      # theta2=fit2$theta
      theta2=fit2$CausalEstimates
    }else{
      theta2=theta1*0
      ind2=1
    }
    vote=normvoting(by=by,bW=bW,theta1=theta1,theta2=theta2,ind1=ind1,ind2=ind2)
    ind1=vote$ind1
    ind2=vote$ind2
    p1=vote$p1
    p2=vote$p2
    error=sqrt(sum(theta11-theta1)^2)
    iter=iter+1
  }
  # cov1=fit1$covtheta
  # cov2=fit2$covtheta
  cov1=fit1$VCovCausalEstimates
  cov2=fit2$VCovCausalEstimates
  
  A=list()
  A$theta1=theta1
  A$theta2=theta2
  A$ind1=ind1
  A$ind2=ind2
  A$delta1=delta1
  A$cov1=cov1
  A$cov2=cov2
  A$p1=p1
  A$p2=p2
  
  if(bic==T){
    fit3=imrp.mrbee.internal(by,bX,Rxx,rxy)
    theta3=fit3$theta
    r1=vec(by[ind1]-bW[ind1,]%*%theta1);r2=vec(by[ind2]-bW[ind2,]%*%theta2);r3=vec(by-bW%*%theta3)
    # bic1=n*log(var(r3))+log(n)*p
    # bic2=n*log(var(r1)^2+var(r2)^2)+log(n)*2*p
    bic1=sum(r3^2)+log(n)*p
    bic2=sum(r1^2*p1)+sum(r2^2*p2)+log(n)*p*2
    if(bic1<bic2){
      A$theta2=A$ind2=A$cov2=0
      A$theta1=theta3
      A$cov1=fit3$covtheta
      A$ind1=c(1:n)
    }
  }
  
  if(aic==T){
    fit3=imrp.mrbee.internal(by,bX,Rxx,rxy)
    theta3=fit3$theta
    r1=vec(by[ind1]-bW[ind1,]%*%theta1);r2=vec(by[ind2]-bW[ind2,]%*%theta2);r3=vec(by-bW%*%theta3)
    bic1=n*log(sd(r3)^2)+2*p
    bic2=log(sd(r1)^2)*length(ind1)+log(sd(r2)^2)*length(ind2)+4*p
    if(bic1<bic2){
      A$theta2=A$ind2=A$cov2=0
      A$theta1=theta3
      A$cov1=fit3$covtheta
      A$ind1=c(1:n)
    }
  }
  wi1=which(names(A)=='theta1')
  wi2=which(names(A)=='cov1')
  wi3=which(names(A)=='delta1')
  names(A)[wi1]='CausalEstimates'
  names(A)[wi2]='VCovCausalEstimates'
  names(A)[wi3]='delta'
  A$method='mix'
  return(A)
}

#' MRBEE with likelihood-weighting of IVs
#'
#' @param pd output of prepData() on the IV set
#' @param max.iter maximum number of iterations to perform until convergence
#' @param max.eps tolerance (smaller values are stricter)
#' @param q tuning parameter in (0,1)
#' @param method method to apply when constructing IV weights
#' @param intercept TRUE/FALSE should an intercept be included in causal estimation  
#' @keywords 
#' @export 
#' @examples
#' MRBEE.MLqe()

MRBEE.MLqe=function(pd,max.iter=15,max.eps=1e-3,q=0.95,method="ordinal",intercept=T){
  by=pd$betaY
  bX=pd$betaX
  Rxx=ssum(pd$UU)/length(pd$betaY)
  rxy=ssum(pd$UV)/length(pd$betaY)
  n=length(by)
  bW=cbind(1,bX)
  Rww=cdiag(0,Rxx)
  rwy=c(0,rxy)
  Hinv=solve(t(bW)%*%bW-n*Rww)
  theta=Hinv%*%(t(bW)%*%by-n*rwy)
  theta[1]=theta[1]*intercept
  theta1=theta*0
  error=sqrt(sum((theta-theta1)^2))
  iter=0
  while(error>max.eps&iter<max.iter){
    theta1=theta
    e=vec(by-bW%*%theta)  
    qweight=mlqeweight(x=e,theta=theta,Rww=Rww,rwy=rwy,q=q,method=method)
    Hinv=solve(t(bW)%*%(bW*qweight)-sum(qweight)*Rww)
    theta=Hinv%*%(t(bW)%*%(by*qweight)-sum(qweight)*rwy)
    theta[1]=theta[1]*intercept
    iter=iter+1
    error=sqrt(sum((theta-theta1)^2))
  }
  E=(-bW*e+kronecker(matrix(1,n,1),t(-vec(Rww%*%theta)+rwy)))*qweight
  V=t(E)%*%E
  covtheta=Hinv%*%V%*%Hinv
  A=list()
  A$CausalEstimates=theta
  A$VCovCausalEstimates=covtheta
  A$qweight=qweight
  A$method='ml'
  return(A)
}

#' A function to perform MR using MRBEE
#'
#' This function performs multivariable/multivariate MR using MRBEE.
#' @param prepDataList the direct output of prepData(), a list of length 5.
#' @param PleioPThreshold The P-value threshold below which a specific IV will be considered horizontally pleiotropic and removed from causal estimation using pleio_test() (see paper).
#' @param FDR logical. If you want to use Benjamini-Hochberg FDR correction instead of defining a P-value threshold for determing horizontal pleiotropy evidence (as in the argument PleioPThreshold), set to TRUE.
#' @param FDR.alpha the uncorrected Type I error rate to be corrected by FDR.
#' @keywords MRBEE
#' @export 
#' @examples
#' MRBEE()

MRBEE=function(prepDataList, PleioPThreshold=0.05, FDR=FALSE, FDR.alpha=0.05) {
  bx=prepDataList$betaX;by=prepDataList$betaY;
  UU=prepDataList$UU;UV=prepDataList$UV;VV=prepDataList$VV
  it=iter.estimator(A=bx,B=by,SigmaUU=UU,SigmaUV=UV,SigmaVV=VV,PleioPThreshold=PleioPThreshold,FDR=FDR,FDR.alpha=FDR.alpha,Outliers=NULL)
  it$method="imrp"
  return(it)
}

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
	Pvalues=1-pchisq(Si, q)
	return(list("Stats"=Si, "Ps"=Pvalues))
}

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

prepData=function(bT,IVInds=1:nrow(bT$EstHarm),oi=1,verbose=TRUE) {
  # bT$R is not ordered in any particular way
  b=bT$EstHarm;s=bT$SEHarm
  betaX=b[IVInds,-oi,drop=FALSE]
  seX=s[IVInds,-oi,drop=FALSE]
  betaY=b[IVInds,oi,drop=FALSE]
  seY=s[IVInds,oi,drop=FALSE]
  inds=1:ncol(b)
  R=bT$R[c(oi,inds[-which(inds==oi)]),c(oi,inds[-which(inds==oi)])]
  # if(any(dim(betaX)!=dim(seX)) | any(dim(betaY)!=dim(seY))) stop("dimensions of betas and standard errors do not all match")
  # if(!all(colnames(betaX)==colnames(seX))) stop("colnames of betaX do not match colnames of seX. Colnames must match and have the same order.")
  # if(!all(colnames(betaY)==colnames(seY))) stop("colnames of betaY do not match colnames of seY. Colnames must match and have the same order.")
  # if(dim(R)[1]!=(ncol(betaX)+ncol(betaY))) stop("dimensions of R do not match number of exposures and outcomes provided")
  # cnR=colnames(R);cnX=colnames(betaX);cnY=colnames(betaY);cnXs=colnames(seX);cnYs=colnames(seY)
  # if(!all(cnR %in% rownames(R))) stop("row and column names of R do not match each other")
  # if(!(all(c(cnX,cnY,cnXs,cnYs) %in% cnR))) stop("please make sure row and column names of matrix R can be found in column names of betaX,betaY,seX,seY")
  # indY=ff(R,cnY); indX=ff(R,cnX); R=R[c(indY,indX),c(indY,indX)]
  m=nrow(betaX);p=ncol(betaX);q=ncol(betaY)
  # outcome was placed in [1,1] position a few lines above
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

#' Plot residuals from MR
#'
#' @param pd output of prepData() on the IV set
#' @param fit fitted model object using one of MRBEE.<x>()
#' @param PleioPThreshold P-value threshold for determining if a specific IV has sufficient evidence of horizontal pleiotropy 
#' @keywords
#' @export
#' @examples
#' residualPlot()

residualPlot=function(pd,fit,PleioPThreshold=0.05) {
  Sp=Spleio(pd,fit)
  bx=pd$betaX;by=c(pd$betaY);byhat=c(cbind(1,bx)%*%fit$CausalEstimates)
  m=length(by)
  softMethods=c('ml')
  hardMethods=c('ipod','mix','imrp','median')
  meth=tolower(fit$method)
  lab1=expression('outcome association ('*hat(bold(alpha))*')')
  lab2=expression('linear prediction of outcome association ('*hat(B)*hat(bold(theta))*')')
  if(meth %in% softMethods) {
    isPleio=fit$qweight
    isPleio=(isPleio-min(isPleio))/(max(isPleio)-min(isPleio)) # rescale to 0-1
    plotdf=data.frame(isPleio, by, byhat, pleioP=Sp$Ps)
    rho=stats::cov.wt(na.omit(cbind(plotdf$by,plotdf$byhat)),wt=fit$qweight,cor=TRUE)$cor[1,2]
    rholab=paste0('linear correlation=',round(rho,2))
    p1=ggplot(plotdf,aes(by,byhat,fill=isPleio,size=-log10(pleioP))) +
      geom_abline(intercept=0,slope=1) +
      # stat_smooth(method='lm',se=FALSE,col='indianred',size=0.6) +
      geom_hline(yintercept=0,color='gray80') +
      geom_vline(xintercept=0,color='gray80') +
      geom_point(alpha=0.5,pch=21) +
      theme_classic() +
      theme(legend.position='bottom') +
      scale_fill_gradient('weight',low='indianred',high='royalblue') +
      # guides(fill=guide_legend(override.aes=list(size=3))) +
      theme(legend.text=element_text(size=9)) +
      labs(x=lab2,y=lab1,subtitle=rholab) +
      guides(size=guide_legend(expression('-log'[10]*'(P'[Pleio]*')')))
  } else if(meth %in% hardMethods) {
    if(meth %in% c('imrp')) fit$delta=c(1:m)*0; fit$delta[fit$PleiotropyIndices]=1
    if(meth=='median') fit$delta=rep(0,m)
    isPleio=fit$delta!=0
    plotdf=data.frame(isPleio=c(isPleio), by, byhat, pleioP=Sp$Ps)
    rho=cor(na.omit(cbind(plotdf$by[!plotdf$isPleio],plotdf$byhat[!plotdf$isPleio])))[1,2]
    rholab=paste0('linear correlation=',round(rho,2))
    p1=ggplot(plotdf,aes(by,byhat,fill=isPleio,size=-log10(pleioP))) +
      geom_abline(intercept=0,slope=1) +
      # stat_smooth(method='lm',se=FALSE,col='indianred',size=0.6,data=plotdf[!isPleio,]) +
      geom_hline(yintercept=0,color='gray80') +
      geom_vline(xintercept=0,color='gray80') +
      geom_point(alpha=0.5,pch=21) +
      theme_classic() +
      theme(legend.position='bottom') +
      scale_fill_manual('',labels=c('Mediation IV','Pleiotropic IV'),
                        values=c('royalblue','gray80')) +
      guides(fill=guide_legend(override.aes=list(size=2))) +
      theme(legend.text=element_text(size=9)) +
      labs(x=lab2,y=lab1,subtitle=rholab) +
      guides(size=guide_legend(expression('-log'[10]*'(P'[Pleio]*')')),
             fill=guide_legend(override.aes=list(size=2,alpha=1)))
  } else {
    stop('did not receive a valid fitted object from MRBEE(), MRBEE.Median(), MRBEE.Mix(), MRBEE.IPOD(), or MRBEE.Median()')
  }
  p1
}

#' A wrapping function to perform an IV-specific test of nonzero horizontal pleiotropy for many IVs at once
#' 
#' @param prepDataList the direct output of prepData(), a list of length 5.
#' @param Outliers vector of indices of IVs to remove
#' @keywords IMRP horizontal pleiotropy test
#' @export 
#' @examples
#' Spleio()

Spleio=function(prepDataList, fit, Outliers=NULL) {
  bx=prepDataList$betaX;by=prepDataList$betaY;
  UU=prepDataList$UU;UV=prepDataList$UV;VV=prepDataList$VV
  p=dim(UU)[1];q=dim(VV)[1]
  SigmaUU0=lapply(array.to.list(UU), function(h) rbind(rep(0, p+1), cbind(rep(0, p), h)))
  SigmaUV0=lapply(array.to.list(UV), function(h) rbind(rep(0, q), matrix(h,nr=p,nc=q)))
  SigmaUU0=simplify2array(SigmaUU0); SigmaUV0=simplify2array(SigmaUV0);
  out=pleio.test(cbind(1,bx),by,fit$CausalEstimates,fit$VCovCausalEstimates,SigmaUU0,SigmaUV0,VV)
  return(out)
}

#' Internal weighting of MRBEE
#'
#' @param pd output of prepData() on the IV set
#' @param weight vector of weights to apply to IVs
#' @param intercept TRUE/FALSE should an intercept be included in causal estimation
#' @keywords
#' @export
#' @examples
#' weight.mrbee.internal()

weight.mrbee.internal=function(pd,weight=rep(1,length(by)),intercept=T){
  by=pd$betaY
  bX=pd$betaX
  Rxx=ssum(pd$UU)/length(pd$betaY)
  rxy=ssum(pd$UV)/length(pd$betaY)
  n=sum(weight)
  bW=cbind(1,bX)
  Rww=cdiag(0,Rxx)
  rwy=c(0,rxy)
  Hinv=solve(t(bW)%*%(bW*weight)-n*Rww)
  theta=vec(Hinv%*%(t(bW)%*%(by*weight)-n*rwy))
  theta[1]=theta[1]*intercept
  e=vec(by-bW%*%theta)
  E=-bW*(e*weight)+kronecker(matrix(1,length(weight),1),t(-vec(Rww%*%theta)+rwy))*weight
  V=t(E)%*%E
  covtheta=Hinv%*%V%*%Hinv
  A=list()
  A$CausalEstimates=theta
  A$VCovCausalEstimates=covtheta
  return(A)
}
