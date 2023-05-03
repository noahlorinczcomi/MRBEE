# source("settingscode.R")
RES1=RES2=matrix(nr=niter,nc=13)
colnames(RES1)=colnames(RES2)=c(
  "MVMR-IVW",#"MRBEE-Mix",
  "MRBEE-IMRP","MRBEE-IMRP-noInt","MVMR-Median","dIVW","MVMR-Lasso","MVMR-Egger",
  "IMRP","MR-Mode","MRMix","MR-CML","MR-Corr","MR-CUE"
)
for(iter in 1:niter) {
  maf=runif(m,0.1,0.9)
  N1=N2=1:n1; N0=(n1-n01+1):(n1-n01+n0); totaln=length(union(N1,N0)) # both exposures always completely overlap
  if(iter%%10==1) {
    G=matrix(0,nr=totaln,nc=m)
    for(snp in 1:m) {
      G[,snp]=rbinom(totaln,2,maf[snp])
      while(sd(G[,snp])==0) G[,snp]=rbinom(totaln,2,maf[snp])
    }
    G=apply(G,2,function(h) (h-mean(h))/sd(h))
  }
  S=matrix(rhox1x2,2,2)+(1-rhox1x2)*diag(2)
  bx12=rmvnorm(m,sigma=S)  
  ss=colSums(bx12^2)
  adjs=h2/ss
  for(i in 1:2) bx12[,i]=sqrt(adjs[i])*bx12[,i]
  # colSums(bx12^2)
  ### exposure model
  U=rnorm(totaln,sd=sqrt((1-h2)*0.15/2))
  epsx1=rnorm(totaln,sd=sqrt(1-h2-var(U)))
  epsx2=rnorm(totaln,sd=sqrt(1-h2-var(U)))
  x1=G%*%bx12[,1]+U+epsx1
  x2=G%*%bx12[,2]+U+epsx2
  # var(G%*%bx)/var(x)
  # outcome model
  epsy=rnorm(totaln,sd=1-sd(x1*theta1+x2*theta2+U))
  uhpix=1:nUHP; chpix=m-c(1:nCHP)+1
  sgns=c();k=0;while(length(sgns)!=m){k=k+1; sgns=c(sgns,ifelse(k%%2==0,1,-1))}
  gammaU=sd(20*bx12%*%c(theta1,theta2))*sgns
  gammaU[-uhpix]=0
  gammaC=-2*bx12%*%c(theta1,theta2)-sd(20*bx12%*%c(theta1,theta2))
  gammaC[-chpix]=0
  if(nUHP==0) gammaU=gammaU*0
  if(nCHP==0) gammaC=gammaC*0
  # plot(bx12%*%c(theta1,theta2),bx12%*%c(theta1,theta2)+gammaU+gammaC)
  #alpha=bx12%*%c(theta1,theta2)+gammaU+gammaC
  y=x1*theta1+G%*%(gammaU+gammaC)+x2*theta2+U+epsy
  y=as.vector(y)
  # for bias-correction
  rhoxy=cor(cbind(x1,x2,y))[3,1:2]
  rhoxx=cor(x1,x2)
  ### separate into non-overlapping (independent) samples
  x2=x2[N2]
  x1=x1[N1]
  y=y[N0]
  G2=G[N2,]
  G1=G[N1,]
  G0=G[N0,]
  # var(x1*theta1+x2*theta2)/var(y) # exposure explains ~12% of variation in outcome
  ### GWAS models
  fit0=biggwas(y,G0); # using G0
  fit1=biggwas(x1,G1); # using G1
  fit2=biggwas(x2,G2); # using G2
  ### MR
  byhat=fit0$est; syhat=fit0$std
  bx1hat=fit1$est; sx1hat=fit1$std
  bx2hat=fit2$est; sx2hat=fit2$std
  bxhat=cbind(bx1hat,bx2hat)
  sxhat=cbind(sx1hat,sx2hat)
  # png(filename="/home/njl96/test.png")
  # plot(bxhat,byhat);abline(a=0,b=theta)
  # points(bxhat[uhpix],byhat[uhpix],col='red')
  # points(bxhat[chpix],byhat[chpix],col='blue')
  # dev.off()
  ### standardization
  # byhat=byhat/syhat
  # bxhat=bxhat/syhat
  uu=m*rhoxx/n1; suu=matrix(c(sum(sx1hat^2),uu,uu,sum(sx2hat^2)),2,2)
  suv=m*rhoxy*c(n01,n02)/sqrt(c(n0*n1,n0*n2))*colMeans(sxhat*syhat)
  svv=m*mean(syhat^2)
  suuArr=array(suu/m,dim=c(2,2,m))
  suvArr=array(suv/m,dim=c(2,1,m))
  svvArr=array(svv/m,dim=c(1,1,m))
  ivw=lm(byhat~bxhat-1,weights=1/syhat^2)
  pd=list(betaX=as.matrix(bxhat),betaY=as.matrix(byhat),UU=suuArr,UV=suvArr,VV=svvArr)
  # mrbmix=MRBEE.Mix(pd,bic=FALSE)
  mrbimrp=MRBEE.IMRP(pd,FDR=TRUE,PleioPThreshold=0.05/m)
  mrbnoint=imrp.mrbee.internal(byhat,bxhat,suu/m,suv/m,max.iter=30,intercept=FALSE); 
  ### MendelianRandomization methods
  ## model fitting
  ################ Only recording causal estimates for the first exposure
  Eivw=coef(ivw)[1]; Sivw=sqrt(diag(vcov(ivw)))[1]
  #Emrbmix=mrbmix$CausalEstimates[2]; Smrbmix=sqrt(mrbmix$VCovCausalEstimates[2,2])
  Emrbimrp=mrbimrp$CausalEstimates[2]; Smrbimrp=sqrt(diag(mrbimrp$VCovCausalEstimates))[2]
  Emrbnoint=mrbnoint$theta[1]; Smrbnoint=sqrt(mrbnoint$covtheta[1,1])
  mvmrobj=MendelianRandomization::mr_mvinput(bx=bxhat,by=byhat,bxse=sxhat,byse=syhat)
  unimrobj=mr_input(bx=bxhat[,1],by=byhat,bxse=sxhat[,1],byse=syhat)
  mrmed=mr_mvmedian(mvmrobj,iterations=1000)
  divw=mr_divw(unimrobj)
  # mrconmix=mr_conmix(unimrobj,CIStep=0.05) # too slow
  # mrmaxlike=mr_maxlik(mrobj) # too slow
  mrlasso=MASS::rlm(bxhat,byhat)
  # mrlasso=mr_lasso(mrobj) # not necessary
  mregg=lm(byhat~bxhat,weights=1/syhat^2)
  imrpdf=data.frame(bxhat=bxhat[,1],byhat=byhat,sxhat=sxhat[,1],syhat=syhat)
  imrp=MR_pleio("byhat","bxhat","syhat","sxhat",imrpdf,0.05/m,rhoxy)
  if(all(is.na(imrp$CausalEstimate))) imrp$CausalEstimate=c(0,1)
  modebased=mr_mbe(unimrobj,iterations=1000)
  mrmix=MRMix(bxhat[,1],byhat,sxhat[,1],syhat,theta_temp_vec=seq(-Emrbimrp-Smrbimrp*5,Emrbimrp+Smrbimrp*5,length.out=20),profile=F)
  ### MVMR-CML is too slow - cannot be ran with>500 IVs
  # Sigma=rbind(suu/m,suv/m); Sigma=cbind(Sigma,c(suv/m,svv/m)); SigmaInv=solve(Sigma)
  # SigmaList=list(); for(si in 1:m) SigmaList[[si]]=SigmaInv
  # MV CML may have some errors and I don't want them to stop progress
  # cml=tryCatch(
  #   MVmr_cML_DP(bxhat,as.matrix(byhat),sxhat,SigmaList,n=min(c(n0,n1,n2)),
  #               K_vec=1:(m/2),num_pert=1,thres=1e-3,maxit=50),
  #   error=function(x) list(BIC_theta=c(NA,NA))
  # )
  # cmlse=tryCatch(
  #   MVcML_SdTheta(bxhat,as.matrix(byhat),SigmaList,cml$BIC_theta,setdiff(1:m,cml$BIC_invalid)),
  #   error=function(x) c(NA,NA)
  # )
  # cml=MVmr_cML_DP(bxhat,as.matrix(byhat),sxhat,SigmaList,n=min(c(n0,n1,n2)),
  #                 K_vec=m/2,num_pert=1,thres=1e-4,maxit=100)
  # cmlse=MVcML_SdTheta(bxhat,as.matrix(byhat),SigmaList,cml$BIC_theta,setdiff(1:m,cml$BIC_invalid))
  nIMRPPleio=min(c(m/2,length(imrp$PleioOutlier)))
  cml=MRcML::mr_cML(b_exp=bxhat[,1],b_out=byhat,se_exp=sxhat[,1],se_out=syhat,n=n0,K_vec=1:nIMRPPleio)
  mrcorropt=list(agm=0.001,bgm=0.001,aal=0.001,bal=0.001,a=1,b=10,maxIter=5000,thin=10,burnin=5000)
  mrcorr=MR.Corr2::MRcorr(gammah=bxhat[,1],Gammah=byhat,se1=sxhat[,1],se2=syhat,opts=mrcorropt)
  mrcueopt=list(agm=0.001,bgm=0.001,atau1=0.001,btau1=0.001,atau2=0.001,btau2=0.001,a=1,b=10,maxIter=5000,thin=10,burnin=5000)
  mrcue=MR.CUE::MRCUEIndep(gammah=bxhat[,1],Gammah=byhat,se1=sxhat[,1],se2=syhat,rho=rhoxy[1],opts=mrcueopt)
  ## extracting estimates
  Emrmed=mrmed@Estimate[1]; Smrmed=mrmed@StdError[1]
  Edivw=divw@Estimate; Sdivw=divw@StdError
  # Emrconmix=mrconmix@Estimate; up=mrconmix@CIUpper; Smrconmix=(up-Emrconmix)/1.96 
  # Emrmaxlike=mrmaxlike@Estimate; Smrmaxlike=mrmaxlike@StdError
  Emrlasso=coef(mrlasso)[1]; Smrlasso=sqrt(diag(vcov(mrlasso)))[1]
  Emregg=coef(mregg)[2]; Smregg=sqrt(diag(vcov(mregg)))[2]
  Eimrp=imrp$CausalEstimate; Simrp=imrp$SdCausalEstimate
  Emodebased=modebased@Estimate; Smodebased=modebased@StdError
  Emrmix=mrmix$theta; Smrmix=mrmix$SE_theta
  # Ecml=cml$BIC_theta[1]; Scml=cmlse[1]
  Ecml=cml$MA_BIC_theta; Scml=cml$MA_BIC_se
  Emrcorr=mean(mrcorr$Beta0res); Smrcorr=sd(mrcorr$Beta0res)
  Emrcue=mrcue$beta.hat; Smrcue=mrcue$beta.se
  ### store results (ignoring MR-CUE)
  RES1[iter,]=c(Eivw,Emrbimrp,Emrbnoint,Emrmed,Edivw,Emrlasso,Emregg,Eimrp,Emodebased,Emrmix,Ecml,Emrcorr,Emrcue)
  RES2[iter,]=c(Sivw,Smrbimrp,Smrbnoint,Smrmed,Sdivw,Smrlasso,Smregg,Simrp,Smodebased,Smrmix,Scml,Smrcorr,Smrcue)
  # save results
  print(iter)
  if(iter%%floor(niter/30)==0) {
    print(colMeans(RES1,na.rm=T))
    # par(mfrow=c(2,1))
    # boxplot(RES1,ylim=c(0,1));abline(h=theta1)
    # boxplot(RES2[,-8]-apply(RES1[,-8],2,function(h) sd(h,na.rm=T)),ylim=c(-1,1));abline(h=0)
    # par(mfrow=c(1,1))
  }
}
fpout=paste0("n",n0,"_pOverlap",p0outcome,"_m",m,"_theta1",theta1,"_theta2",theta2,"_pUHP",propUHP,"_pCHP",propCHP,"FDR005.Rds")
setwd(outputDir)
saveRDS(list(RES1=RES1,RES2=RES2,params=params),file=fpout)
