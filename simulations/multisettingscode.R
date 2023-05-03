library(mvtnorm);library(dplyr);library(ggplot2); library(MASS); library(FDRestimation)
library(MendelianRandomization)
library(MRBEE)
library(MRMix)
library(MR.Corr2)
library(MR.CUE)
library(MRlap)
library(MVMRcML)
library(MR.CUE)
library(IMRP)

n0=n1=n2=20000;
p0outcome=0
p0exposures=1
m=200;
h2=0.05
b=1
theta1=0.5; theta2=0.5
propUHP=0.1
propCHP=0.1

biggwas=function(x,G){
  x=as.vector(x)
  ux=mean(x)
  vx=var(x);vx=as.numeric(vx)
  ug=colMeans(G)
  G=t(t(G)-ug)
  vg=colSums(G^2)
  b=(t(G)%*%(x-ux))/vg
  sdb=(vx-b^2*vg/length(x))/length(x)
  A=list()
  A$est=as.vector(b)
  A$std=as.vector(sqrt(sdb))
  return(A)
}

n01=p0outcome*min(c(n0,n1));n02=p0outcome*min(c(n0,n2));
n12=p0exposures*min(c(n1,n2))
rhox1y=rhox2y=0.25;rhox1x2=0.5
nUHP=floor(propUHP*m)
nCHP=floor(propCHP*m)
niter=500
params=c(n0=n0,n1=n1,n2=n2,m=m,h2=h2,a=a,b=b,theta1=theta1,theta2=theta2,
         propUHP=propUHP,propCHP=propCHP,
         n01=n01,n02=n02,rhox1y=rhox1y,rhox2y=rhox2y,nUHP=nUHP,nCHP=nCHP,niter=niter)
outputDir="./" # change to any existing directory you want
# source("multirunningcode.R")