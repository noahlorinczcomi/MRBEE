
library(MendelianRandomization)
library(MRBEE)
library(FDRestimation)
library(MASS)
ms=c(50,100,250,500,1000)
ps=seq(2,10,1)
niter=5
times=array(dim=c(length(ms),length(ps),7))
for(i in 1:length(ms)) {
  for(j in 1:length(ps)) {
    RES=matrix(nr=niter,nc=7)
    for(iter in 1:niter) {
      bx=rnorm(ms[i]*ps[j],0,1/sqrt(50000));bx=matrix(bx,ms[i],ps[j])
      by=bx%*%rep(0.1,ps[j])+rnorm(ms[i],0,1/sqrt(50000))
      Rho=1/toeplitz(1:(ps[j]+1))
      bxse=sqrt(rchisq(ms[i]*ps[j],50000-1)/(50000-1)^2)
      bxse=matrix(bxse,ms[i],ps[j])
      byse=sqrt(rchisq(ms[i],50000-1)/(50000-1)^2)
      mro=mr_mvinput(bx=bx,by=c(by),bxse=bxse,byse=byse)
      t0=Sys.time()
      mrivw=mr_mvivw(mro)
      res1=as.numeric(Sys.time()-t0)
      # MR-Egger
      t0=Sys.time()
      mregg=mr_mvegger(mro)
      res2=as.numeric(Sys.time()-t0)
      # MR-Lasso
      t0=Sys.time()
      mrlasso=rlm(bx,by)
      res3=as.numeric(Sys.time()-t0)
      # MR-Median
      t0=Sys.time()
      mrmed=mr_mvmedian(mro,iterations=1)
      res4=as.numeric(Sys.time()-t0)
      # MRBEE
      bT=list(R=Rho,Ncor=50000,EstHarm=cbind(by,bx),SEHarm=cbind(byse,bxse))
      pD=prepData(bT)
      t0=Sys.time()
      mrbee=MRBEE.IMRP(pD,FDR=TRUE)
      res5=as.numeric(Sys.time()-t0)
      # MVMR-cML with BIC
      t0=Sys.time()
      cml=tryCatch(mr_mvcML(mro,rho_mat=Rho,K_vec=0:(ms[i]/2),DP=F,n=50000),error=function(x)NA)
      if(is.na(cml)) res6=NA else res6=as.numeric(Sys.time()-t0)
      # MVMR-cML with DP
      t0=Sys.time()
      cml=tryCatch(mr_mvcML(mro,rho_mat=Rho,K_vec=0:(ms[i]/2),DP=T,n=50000),error=function(x)NA)
      if(is.na(cml)) res7=NA else res7=as.numeric(Sys.time()-t0)
      ### store results
      RES[iter,1]=res1
      RES[iter,2]=res2
      RES[iter,3]=res3
      RES[iter,4]=res4
      RES[iter,5]=res5
      RES[iter,6]=res6
      RES[iter,7]=res7
    }
    times[i,j,]=colMeans(RES,na.rm=T)
    print(times)
  }
}
times0=times
saveRDS('times.Rds') # never completed bc MVMR-cML with DP is too slow. only completely finished m through 250 IVs
times=readRDS('C:/Users/njl96/Downloads/times.Rds')
matplot(times[,,7])
times[2:5,,7]=times[2:5,,7]*60 # should be in seconds
times[2,1,7]=times[2,1,7]/60 # should be in minutes
times[3,,1:4]=times[3,,1:4]*60 # should be in seconds

extrap=function(y,x,xex) {
  fit=lm(y~x+I(x^2))
  xp=cbind(1,c(x,xex),c(x,xex)^2)
  cf=ifelse(coef(fit)<1e-6,0,coef(fit))
  xp%*%cf
}
# remove last two rows bc they never finished
times2=array(dim=dim(times)-c(2,0,0))
for(i in 1:dim(times)[3]) times2[,,i]=times[1:3,,i]
times=times2
ms=c(50,100,250,500,1000)[1:3]
ps=seq(2,10,1)
mex=seq(500,5000,length.out=10)
# p=2
r1=matrix(nr=length(c(ms,mex)),nc=dim(times)[3])
for(i in 1:ncol(r1)) r1[,i]=extrap(times[,1,i],ms,mex)
r1[r1<0]=min(r1[r1>0])
# p=10
r2=matrix(nr=length(c(ms,mex)),nc=dim(times)[3])
for(i in 1:ncol(r2)) r2[,i]=extrap(times[,length(ps),i],ms,mex)
# nicer plots
library(ggplot2); library(ggpubr)
methods=c('IVW','MR-Egger','MR-Lasso','MR-Median','MRBEE','MVMR-cML-BIC','MVMR-cML-DP')
r2[5:nrow(r2),2:5]=r2[5:nrow(r2),2:5]*60 # should be in seconds
xb=c(ms,mex)
p1=data.frame(method=rep(methods,each=nrow(r1)),time=c(r1),m=rep(c(ms,mex),ncol(r1))) %>%
  ggplot(aes(m,log(time/60),color=method)) +
  #geom_point() +
  geom_rect(aes(xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf),fill='gray95',alpha=0.7,inherit.aes=FALSE) +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=2,linetype='dashed') +
  geom_hline(yintercept=log(60*5)+0.5,linetype='dashed') +
  geom_text(aes(label=method),size=2,angle=45) +
  #scale_x_continuous(breaks=xb[-2],labels=xb[-2]) +
  theme_classic() +
  theme(legend.position='none') +
  labs(x='number of instrumental variables',
       title='2 exposures',
       y='log time in minutes') +
  lims(x=c(-500,5500),y=c(-10,10)) +
  annotate('text',x=-400,y=-0.5,label='<1 minute',size=3) +
  annotate('text',x=-400,y=2.5,label='7.5 minutes',size=3) +
  annotate('text',x=-400,y=log(60*5)+1,label='5 hours',size=3)
p1
# have to re-do bc was extrapolating to negative values
ts=c()
for(i in 1:length(xb)) {
  mts=c()
  for(iter in 1:500) {
    x=matrix(rnorm(xb[i]*10),xb[i],10)
    y=x%*%rep(0.15,10)+rnorm(xb[i])
    t0=Sys.time()
    fit=lm(y~x)
    mts[iter]=as.numeric(Sys.time()-t0)
  }
  ts[i]=mean(mts)
}
r2[,1]=ts
p2=data.frame(method=rep(methods,each=nrow(r2)),time=c(r2),m=rep(c(ms,mex),ncol(r2))) %>%
  ggplot(aes(m,log(time/60),color=method)) +
  #geom_point() +
  geom_rect(aes(xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf),fill='gray95',alpha=0.7,inherit.aes=FALSE) +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=2,linetype='dashed') +
  geom_hline(yintercept=log(60*5)+0.5,linetype='dashed') +
  geom_text(aes(label=method),size=2,angle=45) +
  #scale_x_continuous(breaks=xb[-2],labels=xb[-2]) +
  theme_classic() +
  theme(legend.position='none') +
  labs(x='number of instrumental variables',
       title='10 exposures',
       y='log time in minutes') +
  lims(x=c(-500,5500),y=c(-12.5,11.5)) +
  annotate('text',x=-400,y=-0.5,label='<1 minute',size=3) +
  annotate('text',x=-400,y=2.5,label='7.5 minutes',size=3) +
  annotate('text',x=-400,y=log(60*5)+1,label='5 hours',size=3)
p2

ggarrange(p1,p2,nrow=1,ncol=2)



















