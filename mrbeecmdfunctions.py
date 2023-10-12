import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas
import numpy
from statsmodels import regression
from numpy import matrix
from scipy import stats
import gzip
import argparse
import os

def usedefault(givenval,newval): # this function exists bc I can't know the default values until I know how many exposures there are
    if givenval=='': # if user wants to use the default, ie let '' be the default value
        ap=newval # make what it should be
    else:
        ap=givenval
    return ap

def flatten_list(_2d_list): # flatten list of lists into a single list
    flat_list=[]
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

def SpleioP(Bhat,Ahat,UU,UV,VV,Thetahat,VThetahat):
    m=len(Bhat[:,0]);q=len(Ahat[0,:])
    I=numpy.eye(q)
    Spleio=[]
    for i in range(0,m):
        K=numpy.kron(I,Bhat[i,:])
        num=Ahat[i,:]-Thetahat.T@Bhat[i,:]
        den=VV+Thetahat.T@UU@Thetahat+K@VThetahat@K.T-2*Thetahat.T@UV
        S=num@numpy.linalg.inv(den)@num.T; S=float([S][0])
        Spleio.append(S)
    P=[1-stats.chi2.cdf(Spleio[x],q) for x in range(0,m)]
    return P

def bootMRBEE(Bhat,Ahat,UU,UV,VV,R,niter=500):
    # R: LD matrix
    m=len(Bhat[:,0]); p=Bhat.shape[1]; ests=numpy.zeros((niter,p))
    for j in range(0,niter):
        ind=random.choices(range(0,m),k=m)
        kBhat=Bhat[ind,:];kAhat=Ahat[ind,:];
        kR1=R[ind,:];kR=kR1[:,ind]; kR=numpy.linalg.pinv(kR) # generalized inverse because not guaranteed to be posdef
        est,V=estMRBEE(kBhat,kAhat,UU,UV,VV,kR)
        estl=[l.tolist() for l in est]
        ests[j,:]=flatten_list(estl)
    cc=pandas.DataFrame(ests).cov()
    ests=numpy.array([numpy.mean(ests[:,x]) for x in range(0,bx.shape[1])]); cc=numpy.array(cc)
    return ests,cc

def estMRBEE(Bhat,Ahat,UU,UV,VV,Rinv,subff=True):
    # Bhat,UU: exposure
    # Ahat,VV: outcome
    # Rinv: inverse of LD matrix
    m=len(Bhat[:,0]);q=len(Ahat[0]);p=Bhat.shape[1]
    est=numpy.linalg.inv(Bhat.T@Rinv@Bhat-m*UU)@(Bhat.T@Rinv@Ahat-m*UV)
    F=(Bhat.T@Rinv@Bhat-m*UU)
    res0=Ahat.squeeze()-(Bhat@est).squeeze()
    ResMat=numpy.zeros((len(res0),len(res0))); numpy.fill_diagonal(ResMat,res0.squeeze())
    if subff==True:
        ff=numpy.zeros((m,p))
        dd=flatten_list([l.tolist() for l in UV-UU@est])
        for i in range(0,ff.shape[0]):
            ff[i,:]=dd
        E=ResMat@Bhat-ff
    else:
        E=ResMat@Bhat
    V=E.T@Rinv@E
    Varest=numpy.linalg.inv(F)@V@numpy.linalg.inv(F)
    return est, Varest

def imrbee(Bhat,Ahat,UU,UV,VV,Rinv,R,PleioPThreshold,boot=False,niter=1000,initial="bee",max_iter=15):
    warnings.filterwarnings('ignore')
    # R: LD matrix
    # Rinv: inverse of LD matrix
    m=len(Bhat[:,0]); p=len(Bhat[0,:])+1; q=len(Ahat[0,:]) # +1 on p bc I will add an intercept without asking
    Bhat_=numpy.column_stack(([1]*Bhat.shape[0],Bhat))
    UU_=numpy.column_stack(([0]*(p-1),UU)); UU_=numpy.row_stack(([0]*p, UU_))
    UV_=numpy.row_stack((0,UV))
    if initial=="robust":
        mod=quantile_regression.QuantReg(Ahat,Bhat_) # includes an intercept
        res=mod.fit(q=0.5)
        est0=res.params.reshape((p,1))
        V0=res.cov_params()
        pleioPs0=SpleioP(Bhat_,Ahat,UU_,UV_,VV,est0,V0)
    else:
        est0,V0=estMRBEE(Bhat_,Ahat,UU_,UV_,VV,Rinv)
        pleioPs0=SpleioP(Bhat_,Ahat,UU_,UV_,VV,est0,V0)
    k=0; thetadiff=1; diffs=[]; kR=Rinv.copy()
    mask=numpy.ones((m,),dtype='bool')
    Outliers0=[i for i in range(0,len(pleioPs0)) if pleioPs0[i]<0.05] # the initial threshold can be stricter than the one used in iteration
    mask[Outliers0]=False
    Outliers0=numpy.array(Outliers0); Outliers=Outliers0.copy().tolist()
    while (k<max_iter) & (thetadiff>(0.0001*(p*q+p))) & (len(Outliers)<(m-10)):   
        k=k+1
        kBhat=Bhat_[mask]
        kAhat=Ahat[mask]
        Rsub=R[mask,:][:,mask]
        kR=numpy.linalg.pinv(Rsub)
        fitkEst,fitkVar=estMRBEE(kBhat,kAhat,UU_,UV_,VV,kR)
        pleioPs0=SpleioP(Bhat_,Ahat,UU_,UV_,VV,fitkEst,fitkVar)
        Outliers=[i for i in range(0,len(pleioPs0)) if pleioPs0[i]<PleioPThreshold]
        mask[Outliers]=0
        if len(Outliers)==kBhat.shape[0]:
            thetadiff=0
            est0=numpy.array([1]*p)*float('nan')
            V0=numpy.zeros((p,p))*float('nan')
            Outliers=list(range(0,m))
            kR=float('nan')
        elif len(Outliers)==0:
            thetadiff=0
            est0,V0,Outliers,kR=fitkEst,fitkVar,Outliers0,kR
        else:
            thetadiff=numpy.sum(abs(fitkEst-est0))
            est0=fitkEst;V0=fitkVar
            diffs.append(thetadiff)
    if (boot is True) & (len(Outliers)<m):
        if len(Outliers)==0:
            kBhat=Bhat.copy(); kAhat=Ahat.copy(); Rsub=R.copy(); kR=numpy.linalg.inv(Rsub)
        else:
            kBhat=numpy.delete(Bhat, (Outliers), axis=0)
            kAhat=numpy.delete(Ahat, (Outliers), axis=0)
            Rsub=subR(R,Outliers)
            kR=numpy.linalg.pinv(Rsub)
        est0,V0=bootMRBEE(kBhat,kAhat,UU,UV,VV,kR,niter=niter)
    warnings.filterwarnings('default')
    return est0, V0, Outliers, kR, k

def genomePleio(BX,BY,BXSE,BYSE,cc,ests,Vests):
    # assumes that BX and BXSE are missing intercept-relevant terms
    BX=numpy.column_stack(([1]*BX.shape[0],BX))
    BXSE=numpy.column_stack(([1]*BX.shape[0],BXSE))
    est=ests.squeeze()
    p=BX.shape[1]; M=BX.shape[0]
    RES=numpy.ones((M,3)) # 3, not p or p+1. ie pleiotropy, exposures, outcome
    progs=numpy.linspace(M/10,M,10); progs=[int(round(progs[_])) for _ in range(0,len(progs))]; progs=numpy.insert(progs,0,0)
    for _ in range(0,M): # I saved a copy of full, non-subsetted bdf, named bdfog, in case user wanted genomewideSpleio
        bxi=BX[_,:].reshape((1,p))
        bxsei=BXSE[_,:].reshape((1,p))
        byi=BY[_].reshape((1,1))
        bysei=BYSE[_].reshape((1,1))        
        D=numpy.diag(numpy.column_stack((bysei,bxsei[:,1:])).squeeze()) # [1:] bc intercept will always be included
        Sigma=D@cc@D
        VV=numpy.array(Sigma.iloc[0,0]).reshape((1,1))
        UU=Sigma.iloc[1:,1:].values
        UV=Sigma.iloc[1:,0].values.reshape((p-1,1))
        # adding intercept-relevant terms
        UU_=numpy.column_stack(([0]*UU.shape[0],UU)); UU_=numpy.row_stack(([0]*UU_.shape[1],UU_))
        UV_=numpy.row_stack((0,UV))
        # association testing
        chistat=((byi-bxi@est)**2)/(VV+est.T@UU_@est+bxi@Vests@bxi.T-2*(est.T@UV_))
        sp=1-stats.chi2.cdf(float(chistat),1)
        exposureP=float(1-stats.chi2.cdf(bxi[:,1:]@numpy.linalg.inv(UU)@bxi[:,1:].T,p)) # joint test for exposures
        outcomeP=float(1-stats.chi2.cdf(byi**2/VV,1)) # test for outcome
        # store results
        RES[_,0]=sp
        RES[_,1]=exposureP
        RES[_,2]=outcomeP
        if _ in progs:
            '{}% of genome-wide Spleio testing complete'.format(round(_/M*100))
    # convert to dataframe
    RESDF=pandas.DataFrame(RES,columns=['SpleioP','JointExposuresP','OutcomeP'])
    return RESDF



