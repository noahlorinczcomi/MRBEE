# source("multisettingscode.R")
RES1=RES2=matrix(nr=niter,nc=13)
colnames(RES1)=colnames(RES2)=c(
  "IVW",#"MRBEE-Mix",
  "MRBEE-IMRP","MRBEE-IMRP-noInt","MR-Median","dIVW","MR-Lasso","MR-Egger",
  "IMRP","MR-Mode","MRMix","MR-CML","MR-Corr","MR-CUE"
)
for(iter in 1:niter) {
  maf=runif(m,0.1,0.9)
  N1=1:n1; N0=(n1-n01+1):(n1-n01+n0); totaln=length(union(N1,N0))
  if(iter%%10==1) {
    G=matrix(0,nr=totaln,nc=m)
    for(snp in 1:m) {
      G[,snp]=rbinom(totaln,2,maf[snp])
      while(sd(G[,snp])==0) G[,snp]=rbinom(totaln,2,maf[snp])
    }
    G=apply(G,2,function(h) (h-mean(h))/sd(h))
  }
  bx=runif(m,-b,b)
  ss=sum(bx^2)
  adjfactor2=h2/ss
  bx=sqrt(adjfactor2)*bx
  # sum(bx^2)
  ### exposure model
  U=rnorm(totaln,sd=sqrt((1-h2)*0.15))
  epsx=rnorm(totaln,sd=sqrt(1-h2-var(U)))
  x=G%*%bx+U+epsx # h2 of exposure
  # var(G%*%bx)/var(x)
  # outcome model
  epsy=rnorm(totaln,sd=1-sd(x*theta+U))
  uhpix=1:nUHP; chpix=m-c(1:nCHP)+1
  sgns=c();k=0;while(length(sgns)!=m){k=k+1; sgns=c(sgns,ifelse(k%%2==0,1,-1))}
  gammaU=sd(20*bx*theta)*sgns
  gammaU[-uhpix]=0
  gammaC=-2*bx*theta-sd(20*bx*theta)
  gammaC[-chpix]=0
  if(nUHP==0) gammaU=gammaU*0
  if(nCHP==0) gammaC=gammaC*0
  epsy=rnorm(totaln,sd=1-sd(x*theta+U))
  alpha=bx*theta+gammaU+gammaC
  # plot(bx,alpha)
  y=G%*%alpha+(U+epsx)*theta+U+epsy
  rhoxy=cor(x,y)
  ### separate into non-overlapping (independent) samples
  x=x[N1]
  y=y[N0]
  G1=G[N1,]
  G0=G[N0,]
  # var(x*theta)/var(y) # exposure explains ~12% of variation in outcome by default
  ### GWAS models
  fit0=biggwas(y,G0); # using G0
  fit1=biggwas(x,G1); # using G1
  ### MR
  byhat=fit0$est; syhat=fit0$std
  bxhat=fit1$est; sxhat=fit1$std
  ### standardization
  # byhat=byhat/syhat
  # bxhat=bxhat/syhat
  suu=sum(sxhat^2); suv=m*rhoxy*n01/sqrt(n0*n1)*mean(sxhat*syhat)
  suuArr=array(suu/m,dim=c(1,1,m))
  suvArr=array(suv/m,dim=c(1,1,m))
  svvArr=array(mean(syhat^2),dim=c(1,1,m))
  ivw=lm(byhat~bxhat-1,weights=1/syhat^2)
  pd=list(betaX=as.matrix(bxhat),betaY=as.matrix(byhat),UU=suuArr,UV=suvArr,VV=svvArr)
  mrbimrp=MRBEE.IMRP(pd,PleioPThreshold=0.05/sqrt(m),FDR=TRUE,FDR.alpha=0.05)
  mrbnoint=imrp.mrbee.internal(byhat,bxhat,suu/m,suv/m,max.iter=30,intercept=FALSE); 
  ### MendelianRandomization methods
  ## model fitting
  Eivw=coef(ivw); Sivw=sqrt(diag(vcov(ivw)))[1]
  Emrbimrp=mrbimrp$CausalEstimates[2]; Smrbimrp=sqrt(diag(mrbimrp$VCovCausalEstimates))[2]
  Emrbnoint=mrbnoint$theta; Smrbnoint=sqrt(mrbnoint$covtheta)
  mrobj=mr_input(bx=bxhat,by=byhat,bxse=sxhat,byse=syhat)
  mrmed=mr_median(mrobj,iterations=500)
  divw=mr_divw(mrobj)
  # mrconmix=mr_conmix(mrobj,CIStep=0.5) # too slow
  # mrmaxlike=mr_maxlik(mrobj) # too slow
  mrlasso=MASS::rlm(bxhat,byhat)
  # mrlasso=mr_lasso(mrobj) # not necessary
  mregg=lm(byhat~bxhat,weights=1/syhat^2)
  imrpdf=data.frame(bxhat=bxhat,byhat=byhat,sxhat=sxhat,syhat=syhat)
  imrp=MR_pleio("byhat","bxhat","syhat","sxhat",imrpdf,0.05/m,rhoxy)
  if(all(is.na(imrp$CausalEstimate))) imrp$CausalEstimate=c(0,1)
  modebased=mr_mbe(mrobj,iterations=500)
  mrmix=MRMix(bxhat,byhat,sxhat,syhat,theta_temp_vec=seq(-Emrbimrp-Smrbimrp*5,Emrbimrp+Smrbimrp*5,length.out=20),profile=F)
  nIMRPPleio=min(c(m/2,length(imrp$PleioOutlier)*2))
  cml=MRcML::mr_cML(b_exp=bxhat,b_out=byhat,se_exp=sxhat,se_out=syhat,n=n0,K_vec=1:nIMRPPleio,maxit=100)
  mrcorropt=list(agm=0.001,bgm=0.001,aal=0.001,bal=0.001,a=1,b=10,maxIter=5000,thin=10,burnin=5000)
  mrcorr=MR.Corr2::MRcorr(gammah=bxhat,Gammah=byhat,se1=sxhat,se2=syhat,opts=mrcorropt)
  mrcueopt=list(agm=0.001,bgm=0.001,atau1=0.001,btau1=0.001,atau2=0.001,btau2=0.001,a=1,b=10,maxIter=5000,thin=10,burnin=5000)
  mrcue=MR.CUE::MRCUEIndep(gammah=bxhat,Gammah=byhat,se1=sxhat,se2=syhat,rho=rhoxy,opts=mrcueopt)
  ## extracting estimates
  Emrmed=mrmed@Estimate; Smrmed=mrmed@StdError
  Edivw=divw@Estimate; Sdivw=divw@StdError
  # Emrconmix=mrconmix@Estimate; up=mrconmix@CIUpper; Smrconmix=(up-Emrconmix)/1.96
  # Emrmaxlike=mrmaxlike@Estimate; Smrmaxlike=mrmaxlike@StdError
  Emrlasso=coef(mrlasso); Smrlasso=sqrt(diag(vcov(mrlasso)))[1]
  Emregg=coef(mregg)[2]; Smregg=sqrt(diag(vcov(mregg)))[2]
  Eimrp=imrp$CausalEstimate; Simrp=imrp$SdCausalEstimate
  Emodebased=modebased@Estimate; Smodebased=modebased@StdError
  Emrmix=mrmix$theta; Smrmix=mrmix$SE_theta
  Ecml=cml$MA_BIC_theta; Scml=cml$MA_BIC_se
  Emrcorr=mean(mrcorr$Beta0res); Smrcorr=sd(mrcorr$Beta0res)
  Emrcue=mrcue$beta.hat; Smrcue=mrcue$beta.se
  ### store results (ignoring MR-CUE)
  RES1[iter,]=c(Eivw,Emrbimrp,Emrbnoint,Emrmed,Edivw,Emrlasso,Emregg,Eimrp,Emodebased,Emrmix,Ecml,Emrcorr,Emrcue)
  RES2[iter,]=c(Sivw,Smrbimrp,Smrbnoint,Smrmed,Sdivw,Smrlasso,Smregg,Simrp,Smodebased,Smrmix,Scml,Smrcorr,Smrcue)
  #if(iter%%10==0) { boxplot(RES);abline(h=theta,col="red"); print(iter) }
  # save results
  print(iter)
  if(iter%%floor(niter/30)==0) {
    print(colMeans(RES1,na.rm=T))
    # boxplot(RES1,ylim=c(0.25,0.75));abline(h=theta)
  }
}
fpout=paste0("n",n0,"_pOverlap",p0,"_m",m,"_theta",theta,"_pUHP",propUHP,"_pCHP",propCHP,"FDR005.Rds")
setwd(outputDir)
saveRDS(list(RES1=RES1,RES2=RES2,params=params),file=fpout)
