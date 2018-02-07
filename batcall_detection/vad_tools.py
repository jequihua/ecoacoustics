#import numpy as np
import numpy.matlib as ml
from numpy import size,maximum,minimum,remainder,sqrt,exp,log,max,min,ones,zeros,round,int,array,nonzero
from numpy import full,inf,NaN,logical_not,logical_and,logical_or,floor,pi, count_nonzero, delete,concatenate

from scipy.special import iv as besseli
from scipy.signal  import spectrogram as spec

#from numpy.lib import stride_tricks
def getvadparam(dT=5e-3,theta=0.5):
    # PARAMETERS FOR ADJUSTED FOR ULTRASONIC BAT CALLS 
    # vad.dt  = dT;             % true frame increment time
    # vad.pr  = theta;            % Speech probability threshold
    # vad.gx  = db2pow(30);     % maximum posterior SNR = 30dB
    # vad.gz  = db2pow(-40);    % minimum posterior SNR = -40dB

    # vad.ne  = 0;              % noise estimation: 0=min statistics, 1=MMSE [0]

    # vad.ge  = 1;                % xi estimation: 0= Itakura Saito (ISD), 1 = Decision Directed (DD) [1]
    # vad.ta  = -dT/log(0.98);    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [2]
    # vad.xn  = 0;                % minimum prior SNR = -Inf dB

    # vad.hmm = 1;                % HMM-Based Hang Over [0]
    # vad.ts  = 0.04;             % mean talkspurt length (100 ms)
    # vad.tn  = 0.005;             % mean silence length (50 ms)


    # % Estimate noise spectrum using minimum statistics
    # ne.taca   =-dT/log(0.95);       % smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11)
    # ne.tamax  =-dT/log(0.96);      % max smoothing time constant in (3) = -tinc/log(0.96)
    # ne.taminh =-dT/log(0.3);       % min smoothing time constant (upper limit) in (3) = -tinc/log(0.3)
    # ne.tpfall =0.064/10;           % time constant for P to fall (12)
    # ne.tbmax  =-dT/log(0.8);             % max smoothing time constant in (20) = -tinc/log(0.8)
    # ne.qeqmin =2;                  % minimum value of Qeq (23)
    # ne.qeqmax =14;                 % max value of Qeq per frame
    # ne.av     =2.12;               % fudge factor for bc calculation (23 + 13 lines)
    # ne.td     =0.01;              % time to take minimum over
    # ne.nu     =16;                  % number of subwindows
    # ne.qith   =[0.03 0.05 0.06 Inf];% noise slope thresholds in dB/s
    # ne.nsmdb  =[47 31.4 15.7 4.1];    

    vad = dict(
    dt = dT, 
    pr = theta,
    gx = 10**(30/10),gz = 10**(-40/10),
    ne = 0,
    ge = 1,
    ta = -dT/log(0.98),
    xn = 0,
    hmm = 1,
    ts = 0.04,
    tn = 0.005)


    ne = dict(
    taca   = -dT/log(0.95),       
    tamax  = -dT/log(0.96),     
    taminh = -dT/log(0.3),       
    tpfall = 0.064/10,           
    tbmax  = -dT/log(0.8),            
    qeqmin = 2.0,                  
    qeqmax = 14.0,                 
    av     = 2.12,               
    td     = 0.01,                 
    nu     = 16.0,                  
    qith   = ml.asarray([0.03, 0.05, 0.06, inf]),
    nsmdb  = ml.asarray([47.0, 31.4, 15.7, 4.1]))
    return(vad,ne)
    
def PSDcompute(x,fs,window='hamming',nperseg=1024,overlap=0.75):
    
    # Audio normalization
    S0    = x.nonzero()[0][0]
    T0    = S0/fs    
    x     = x[S0:]/x.max()*0.9
     
    # SFTF calculation        
    noverlap = floor(overlap*nperseg)
    F,T,P    = spec(x,fs,window,nperseg,noverlap,nfft=None,
                    detrend='constant', return_onesided=True, scaling='density', 
                    axis=-1, mode='psd')
    # STFT filtering
    iF = logical_and(12.5e3<F,F<300e3)
    f  = F[iF]
    t  = T+T0
    p  = P[iF,:]

    dt = T[2]-T[1]
    df = F[2]-F[1]
    return(p,f,t,dt,df)
    
def vadcompute(X,N,vad,ne):  
    
    dt   = vad['dt']          
    ta   = vad['ta']
    gmax = vad['gx']          
    xn   = vad['xn']          
    gz   = vad['gz']          
    tn   = vad['tn']
    ts   = vad['ts']
    theta= vad['pr']
    hmm  = vad['hmm']

    a  = exp(-dt/ta)        # SNR smoothing coefficient
    kk = sqrt(2*pi)      # sqrt(8)*Gamma(1.5) - required constant    
    
    a01 = dt/tn                 # a01=P(silence->speech)
    a00 = 1-a01                 # a00=P(silence->silence)
    a10 = dt/ts                 # a10=P(speech->silence)
    a11 = 1-a10                 # a11=P(speech->speech)
    
    P01 = a10/a01               # P(silence)/P(speech)
    logP01= log(P01)

    (I,J)   = X.shape
    #Pt      = X.transpose()    
    #if vad['ne']==1:
    #    dp = estnoiseg(Pt,dt,ne) # estimate the noise using MMSE
    #elif vad['ne']==0:         
    #    dp = estnoisem(Pt,dt,ne) # estimate the noise using minimum statistics

    if vad['ge']==1:    # Decision Directed
        xu      = 1     # dummy unsmoothed SNR from previous frame [2](53)++
    elif vad['ge']==0:  # Maximum Likelihood        
        pass
    gamma   = maximum(minimum(X/N,gmax),gz)
    logl    = ml.asarray(zeros((I,J)))
    logL    = ml.asarray(zeros((1,J)))
    PR      = logL.copy()                # create space for prob ratio vector
    LR      = logL.copy()
    
    
    for j in range(J):
        gammak=gamma[:,j];

        #  Decision Directed
        xik = a*xu+(1-a)*maximum(gammak-1,xn)    # xi = smoothed a priori SNR [2](53)
        v  = 0.5*xik*gammak/(1+xik)
        gi =(0.277+2*v)/gammak                 # accurate to 0.02 dB for v>0.5
        mv = v<0.5
        if ml.any(mv):                               # only calculate Bessel functions for v<0.5
            vmv   = v[mv]
            eq1     = (0.5 + vmv)*besseli(0,vmv) + vmv*besseli(1,vmv)
            eq2     = gammak[mv]*exp(vmv)
            gi[mv]  = kk*sqrt(vmv)*eq1/eq2  # [2](7)
        
        xu = gammak*(gi ** 2)                     # unsmoothed prior SNR % [2](53)
        #   Maximum Likelihood
        #   xu=gammak.*gi.^2;
        
        loglj    = (xik*gammak/(1+xik)) - log(1+xik)
        logl[:,j]= loglj.copy()
        logL[0,j]  = loglj.sum()/I

        if hmm==1:
            
            if j==0:
                Gamma=- logP01 + logL[0,j]            
            else:                
                Gamma = (a01+a11*Gamma)/(a00+a10*Gamma)*exp(logL[0,j])
            logLn = logP01 + log(Gamma)
            LR[0,j] = logLn
            PR[0,j] = 2.0 / (1 + exp(- logLn)) - 1
        elif hmm==0:
            LR[0,j] = logL[0,j]
            PR[0,j] = 2.0 / (1 +exp(- logL[j])) - 1

    return LR,logl

def peakdetection(S,Sxx,theta,t,f,vad,ne):
    Signal = S[0]
    (I,J)   = Sxx.shape
    # Peak detection
    D      = (Signal >= theta)*1 
    Dif    = ml.diff(D)
    onset         = full(J,False)
    offset        = full(J,False)
    onset [1:J]   = Dif== 1
    offset[0:J-1] = Dif==-1
    
    if  count_nonzero(onset)  < count_nonzero(offset):
        onset[0]  = True
    elif count_nonzero(onset) > count_nonzero(offset):
        offset[-1]= True
    
    Onset = nonzero(onset)[0]
    Offset= nonzero(offset)[0]


    if vad['ne'] == 0:          #delete intial false detection
        if ml.any(onset):
            if Onset[0] < round(ne['td'] / vad['dt']):
                D[Onset[0]:Offset[0]] = False
                Onset  = delete(Onset,0)
                Offset = delete(Offset,0)
   
    
    P     = Onset.size
    if P>0:
        iPmax = zeros(P)
        Pmax  = full ([P],NaN)
        
        for p in range(P):
            if Onset[p] == Offset[p] :
                signal      = Signal[Onset[p]]
            else :
                signal      = Signal[Onset[p]:Offset[p]]
            Pmax[p]     = signal.max()
            iPmax[p]    = signal.argmax() + Onset[p]

        iPx   = array([Onset,iPmax,Offset])
        Px    = Pmax
    
        SxxPy = Sxx[:,iPmax.astype('int')]
        Py     = SxxPy.max(0)
        iPy    = SxxPy.argmax(0)
        
        iPt = t[iPx.astype('int')].T
        iPf = f[iPy.astype('int')].T
        
        out = dict(
                    NumberOfDetections    = P,
                    DetectionTimeInt      = array([iPt[:,0], iPt[:,2]]),
                    DetectionTimeLength   = iPt[:,2] - iPt[:,0],
                    DetectionTimePoint    = iPt[:,1],
        
                    DetectionFreqInt      = Py,
                    DetectionFreqPoint    = iPf,
        
                    DetectionSmpPoint     = array([iPx.T[:,1], iPy]),
                    DetectionPowPoint     = array([Px,Py]),
                    DetectionSignal       = array([D])
                    )
                    

    return(out)
    
    
def estnoisem(yf=None,dt=None,ne=None):
    #ESTNOISEM - estimate noise spectrum using minimum statistics
    # Usage:    ninc=round(0.016*fs);   # frame increment [fs=sample frequency]
    #           ovf=2;                  # overlap factor
    #           f=rfft(enframe(s,hanning(ovf*ninc,'periodic'),ninc),ovf*ninc,2);
    #           f=f.*conj(f);           # convert to power spectrum
    #           x=estnoisem(f,ninc/fs); # estimate the noise power spectrum
    
    # Inputs:
    #   yf      input power spectra (one row per frame)
    #   tz      frame increment in seconds
    #           Alternatively, the input state from a previous call (see below)
    #   pp      algorithm parameters [optional]
    
    # Outputs:
    #   x       estimated noise power spectra (one row per frame)
    #   zo      output state
    #   xs      estimated std error of x (one row per frame)
    #           xs seems often to be an underestimate by a factor of 2 or 3
    
    # The algorithm parameters are defined in reference [1] from which equation
    # numbers are given in parentheses. They are as follows:
    
    #        pp.taca      # (11): smoothing time constant for alpha_c [0.0449 seconds]
    #        pp.tamax     # (3): max smoothing time constant [0.392 seconds]
    #        pp.taminh    # (3): min smoothing time constant (upper limit) [0.0133 seconds]
    #        pp.tpfall    # (12): time constant for P to fall [0.064 seconds]
    #        pp.tbmax     # (20): max smoothing time constant [0.0717 seconds]
    #        pp.qeqmin    # (23): minimum value of Qeq [2]
    #        pp.qeqmax    # max value of Qeq per frame [14]
    #        pp.av        # (23)+13 lines: fudge factor for bc calculation  [2.12]
    #        pp.td        # time to take minimum over [1.536 seconds]
    #        pp.nu        # number of subwindows to use [3]
    #        pp.qith      # Q-inverse thresholds to select maximum noise slope [0.03 0.05 0.06 Inf ]
    #        pp.nsmdb     # corresponding noise slope thresholds in dB/second   [47 31.4 15.7 4.1]
    
    # Example use:      y=enframe(s,w,ni);                  # divide speech signal s(n) into
    #                                                       # overlapping frames using window w(n)
    #                   yf=rfft(y,nf,2);                    # take fourier transform
    #                   dp=estnoisem(yf.*conj(yf),tinc);    # estimate the noise
    
    # If convenient, you can call estnoisem in chunks of arbitrary size. Thus the following are equivalent:
    
    #                   (a) dp=estnoisem(yp(1:300),tinc);
    
    #                   (b) [dp(1:100),z]=estnoisem(yp(1:100),tinc);
    #                       [dp(101:200),z]=estnoisem(yp(101:200),z);
    #                       [dp(201:300),z]=estnoisem(yp(201:300),z);
    
    # This is intended to be a precise implementation of [1] with Table III
    # replaced by the updated table 5 from [2]. The only deliberate algorithm
    # change is the introduction of a minimum value for 1/Qeq in equation (23).
    # This change only affects the first few frames and improves the
    # convergence of the algorithm. A minor improveemnt was reported in [3] but
    # this has not yet been included.
    
    # Refs:
    #    [1] Rainer Martin.
    #        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
    #        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.
    #    [2] Rainer Martin.
    #        Bias compensation methods for minimum statistics noise power spectral density estimation
    #        Signal Processing, 2006, 86, 1215-1229
    #    [3] Dirk Mauler and Rainer Martin
    #        Noise power spectral density estimation on highly correlated data
    #        Proc IWAENC, 2006
    
    #	   Copyright (C) Mike Brookes 2008
    #      Version: $Id: estnoisem.m 1718 2012-03-31 16:40:41Z dmb $
    
    #   VOICEBOX is a MATLAB toolbox for speech processing.
    #   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
    
    #################################################################################
    #   This program is free software; you can redistribute it and/or modify
    #   it under the terms of the GNU General Public License as published by
    #   the Free Software Foundation; either version 2 of the License, or
    #   (at your option) any later version.
    
    #   This program is distributed in the hope that it will be useful,
    #   but WITHOUT ANY WARRANTY; without even the implied warranty of
    #   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #   GNU General Public License for more details.
    
    #   You can obtain a copy of the GNU General Public License from
    #   http://www.gnu.org/copyleft/gpl.html or by writing to
    #   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
    #################################################################################
    yf    = yf.T
    nr,nrf= yf.shape
    N     = zeros((nr,nrf))
    #xs    =zeros(nr,nrf)
    tinc  = dt
    nrcum = 0.
    
    # default algorithm constants
    
    # unpack parameter structure
    taca    = ne['taca']        # smoothing time constant for alpha_c = -tinc/log(0['7) in equ (11)
    tamax   = ne['tamax']       # max smoothing time constant in (3) = -tinc/log(0['96)
    taminh  = ne['taminh']      # min smoothing time constant (unppr limit) in (3) = -tinc/log(0['3)
    tpfall  = ne['tpfall']      # time constant for P to fall (12)
    tbmax   = ne['tbmax']       # max smoothing time constant in (20) = -tinc/log(0['8)
    qeqmin  = ne['qeqmin']      # minimum value of Qeq (23)
    qeqmax  = ne['qeqmax']      # max value of Qeq per frame
    av      = ne['av']          # fudge factor for bc calculation (23 + 13 lines)
    td      = ne['td']          # time to take minimum over
    nu      = ne['nu']          # number of subwindows
    qith    = ne['qith']        #noise slope thresholds in dB/s
    nsmdb   = ne['nsmdb']       # maximum permitted +ve noise slope in dB/s    

    
    # derived algorithm constants
    
    aca    =exp(- tinc/taca)     
    acmax  =aca                  
    amax   =exp(- tinc/tamax)    
    aminh  =exp(- tinc/taminh)    
    bmax   =exp(- tinc/tbmax)
    snrexp =-tinc/tpfall
    nv     = round(td/(tinc*nu))    
    if nv < 4:
        nv =4.
        nu = maximum(round(td/(tinc*nv)),1)
        nu = int(nu)   
    nd        = nu*nv    
    md,hd,xx     = mhvals(nd)   
    mv,hv,xx     = mhvals(nv)
    nsms      = 10.0 ** (nsmdb*nv*tinc/10.)
    qeqimax   =1/qeqmin
    qeqimin   =1/qeqmax
    
    if logical_not(nrcum):
        P        = array([yf[1,:]])
        ac       = 1.
        sn2      = P
        pb       = P
        pb2      = pb ** 2.
        pminu    = P
        actmin   = full([1,nrf],inf)
        actminsub= actmin
        subwc    = nv
        actbuf   = full([nu,nrf],inf)
        ibuf     = 0
        lminflag = zeros((1,nrf))
    # loop for each frame    
    for t in range(nr):
        yft     = array([yf[t,:]])        # noise speech power spectrum
        acb     = (1+(P.sum()/yft.sum()-1)**2)**(-1)  # alpha_c-bar(t)  (9)
        ac      = aca*ac+(1-aca)*maximum(acb,acmax)       # alpha_c(t)  (10)
        ah      = amax*ac*(1+(P/sn2-1)**2)**(-1)    # alpha_hat: smoothing factor per frequency (11)
        snr     = P.sum()/sn2.sum()
        ah      = maximum(ah,minimum(aminh,snr**snrexp))       # lower limit for alpha_hat (12)

        P       = ah*P+(1-ah)*yft            # smoothed noisy speech power (3)
        b       = minimum(ah**2,bmax)              # smoothing constant for estimating periodogram variance (22 + 2 lines)
        pb      = b*pb + (1-b)*P            # smoothed periodogram (20)
        pb2     = b*pb2 + (1-b)*(P**2)       # smoothed periodogram squared (21)

        qeqi    = maximum(minimum((pb2-pb**2)/(2*sn2**2),qeqimax),qeqimin/((t+1.)+nrcum))   # Qeq inverse (23)
        qiav    = qeqi.sum()/nrf             # Average over all frequencies (23+12 lines) (ignore non-duplication of DC and nyquist terms)
        bc      = 1+av*sqrt(qiav)             # bias correction factor (23+11 lines)
        bmind   = 1+2*(nd-1)*(1-md)/(qeqi**(-1)-2*md)      # we use the simplified form (17) instead of (15)
        bminv   = 1+2*(nv-1)*(1-mv)/(qeqi**(-1)-2*mv)      # same expression but for sub windows
        kmod    = bc*P*bmind < actmin        # Frequency mask for new minimum
        if kmod.any():
                actmin[kmod]   = bc*P[kmod]*bmind[kmod]
                actminsub[kmod]= bc*P[kmod]*bminv[kmod]
        
        if   subwc > 1 and subwc < nv:
            lminflag = logical_or(lminflag,kmod)
            pminu    = minimum(actminsub,pminu)
            sn2      = pminu
        elif subwc>=nv:                         # end of buffer - do a buffer switch
                    ibuf          = 0+remainder(ibuf,nu)        # increment actbuf storage pointer
                    actbuf[ibuf,:]= actmin      # save sub-window minimum
                    pminu         = actbuf.min(0)

                    i     = nonzero(qiav<qith)[0]
                    nsm   = nsms[i[0]]             # noise slope max
                    nkmod = logical_not(kmod)
                    lmin  = logical_and(logical_and(logical_and(lminflag,nkmod),actminsub<nsm*pminu), actminsub>pminu)
                    if lmin.any():
                        pminu [lmin]   =actminsub[lmin]
                        actbuf[:,lmin] =ones(nu)*pminu[lmin]
                    lminflag[:] = 0.
                    actmin[:]   = inf
                    subwc       = 0.
        subwc=subwc + 1.
        N[t,:]= sn2
        #qisq  = sqrt(qeqi)
        # empirical formula for standard error based on Fig 15 of [2]
        #xs(t,:)=sn2.*sqrt(0.266*(nd+100*qisq).*qisq/(1+0.005*nd+6/nd)./(0.5*qeqi.^(-1)+nd-1))
        
    return N.T
    
def mhvals(D=None):
    # Values are taken from Table 5 in [2]
    #[2] R. Martin,"Bias compensation methods for minimum statistics noise power
    #               spectral density estimation", Signal Processing Vol 86, pp1215-1229, 2006.
    # approx: plot(d.^(-0.5),[m 1-d.^(-0.5)],'x-'), plot(d.^0.5,h,'x-')
    
    DMH=      array([[1,0,0],
                     [2,0.26,0.15],
                     [5,0.48,0.48],
                     [8,0.58,0.78],
                     [10,0.61,0.98],
                     [15,0.668,1.55],
                     [20,0.705,2],
                     [30,0.762,2.3],
                     [40,0.8,2.52],
                     [60,0.841,3.1],
                     [80,0.865,3.38],
                     [120,0.89,4.15],
                     [140,0.9,4.35],
                     [160,0.91,4.25],
                     [180,0.92,3.9],
                     [220,0.93,4.1],
                     [260,0.935,4.7],
                     [300,0.94,5]])

    i = nonzero(D <= DMH[:,0])[0]
    if size(i)== 0:
        i = int(DMH.shape[0])
        j = i
    else:
        i=i[0]
        j=i - 1
    if D == DMH[i,0]:
        M=DMH[i,1]
        H=DMH[i,2]
    else:
        qj=sqrt(DMH[i - 1,0])
        qi=sqrt(DMH[i,0])
        qq=sqrt(D)
        H=DMH[i,2] + (qq - qi)*(DMH[j,2] - DMH[i,2]) / (qj - qi)
        M=DMH[i,1] + (qi*qj/qq - qj)*(DMH[j,1] - DMH[i,1]) / (qi - qj)
    return M,H,D

