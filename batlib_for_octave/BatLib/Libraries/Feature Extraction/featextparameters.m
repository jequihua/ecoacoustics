function FE = featextparameters(inFE,Data)

    det.flim        = Data.Ffilter;
    det.wn          = Data.WindowSize;
    det.wt          = det.wn./Data.Fs;
    det.S           = Data.SampLim;
    det.Smax        = Data.TotalSamples;
    det.dS          = abs(diff(Data.SampLim));
    det.Tp          = Data.Tpeak;
    det.Sp          = det.Tp.*Data.Fs;
    det.dSp         = det.wn;
   
   
    stft.fplim    = Data.Ffilter;
    stft.wt       = inFE.WindowTime;
    stft.ovp      = inFE.Overlap;
    stft.wtype    = inFE.WindowType;
    eval(['stft.whandle = @' stft.wtype ';'])
    stft.df       = inFE.deltaF;
    stft.dt       = inFE.deltaT;
    stft.aframe   = stft.df*stft.dt;
    stft.amax     = inFE.MaximumArea/stft.aframe;
    stft.amin     = inFE.MinimumArea/stft.aframe;    
    stft.thetha   = inFE.Threshold;
    stft.delta    = inFE.AdjustmentStep;
    stft.mode     = inFE.FollowFreqPeak;
    
    for n = 1: length(Data.Fs)
    stft.wn(n,1)      = nearest2n (stft.wt,Data.Fs(n));
    end
    stft.ov     = round (stft.wn*stft.ovp);    
    stft.dS     = round(stft.dt.*Data.Fs);
    
    Sp = det.Sp;
    Si = det.S(:,1);
    Sf = det.S(:,2);

    S0i = round(Sp - max(Sp-Si,stft.dS/2) - det.dSp/2);
    S0i = max(S0i,1);
    S0f = round(Sp + max(Sf-Sp,stft.dS/2) + det.dSp/2);
    S0f = min(S0f,det.Smax);
    det.S0 = [S0i S0f];
    det.T0 = det.S0./[Data.Fs Data.Fs];
    det.s  = [Si Sf] - [S0i S0i] +1;

    stft.sp     = Sp - S0i +1;
    stft.tp     = stft.sp./Data.Fs;
    si          = stft.tp - stft.wt/2 - det.wt/2;
    sf          = stft.tp + stft.wt/2 + det.wt/2;
    stft.tplim  = [si sf];
    
    FE.stft = stft;
    FE.det  = det; 
    
    