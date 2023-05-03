library(mvtnorm);library(dplyr);library(ggplot2); library(MASS); library(FDRestimation)
library(MendelianRandomization)
library(MRBEE)
library(MRMix)
library(MR.Corr2)
library(MR.CUE)
library(MRcML)
library(IMRP)
library(MR.CUE)

n0=n1=20000;
p0=1
m=100;h2=0.05 # m=100, 500
a=0;b=1 # range of uniform dist to draw exopsure effect sizes from
theta=0.5
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

n01=p0*min(c(n0,n1));
nUHP=floor(propUHP*m)
nCHP=floor(propCHP*m)
niter=500
params=c(n0=n0,n1=n1,m=m,h2=h2,a=a,b=b,theta=theta,
         propUHP=propUHP,propCHP=propCHP,
         n01=n01,nUHP=nUHP,nCHP=nCHP,niter=niter)
outputDir="./" # change to any existing directory you want
# source("unirunningcode.R")
