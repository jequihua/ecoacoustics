function FEData = featextparameters(inFE,Data)

    FEData.det.flim        = Data.Ffilter;
    FEData.det.wn          = Data.WindowSize;
    for n = 1: Data.N
    FEData.det.wt(n)      = nearest2n (Data.WindowSize(n),Data.Fs(n));
    end
    FEData.det.S           = Data.SampLim;
    FEData.det.Smax        = Data.TotalSamples;
    FEData.det.dS          = abs(diff(Data.SampLim));
    FEData.det.Tp          = Data.Tpeak;
    FEData.det.Sp          = det.Tp.*Data.Fs;
    FEData.det.dSp         = det.wn;
   
   
    FEData.stft.fplim    = Data.Ffilter;
    FEData.stft.wt       = inFE.WindowTime;
    FEData.stft.ovp      = inFE.Overlap;
    FEData.stft.wtype    = inFE.WindowType;
    eval(['frame.stft.whandle = @' stft.wtype ';'])
    FEData.stft.df       = inFE.deltaF;
    FEData.stft.dt       = inFE.deltaT;
    FEData.stft.aframe   = stft.df*stft.dt;
    FEData.stft.amax     = inFE.MaximumArea/stft.aframe;
    FEData.stft.amin     = inFE.MinimumArea/stft.aframe;    
    FEData.stft.thetha   = inFE.Threshold;
    FEData.stft.delta    = inFE.AdjustmentStep;
    FEData.stft.mode     = inFE.FollowFreqPeak;
    
    for n = 1: Data.N
    FEData.stft.wn(n)      = nearest2n (stft.wt,Data.Fs(n));
    end
    FEData.stft.ov     = round (wn*stft.ovp);    
    FEData.stft.dS     = round(stft.dt.*Data.Fs);
    
    Sp = det.Sp;
    Si = det.S(:,1);
    Sf = det.S(:,2);

    S0i = round(Sp - max(Sp-Si,stft.dS/2) - det.dSp/2);
    S0i = max(S0i,1);
    S0f = round(Sp + max(Sf-Sp,stft.dS/2) + det.dSp/2);
    S0f = min(S0f,det.Smax);
    FEData.det.S0 = [S0i S0f];
    FEData.det.s  = [Si Sf] - [S0i S0i] +1;

    FEData.stft.sp     = Sp - S0i +1;
    FEData.stft.tp     = FEData.stft.sp/Data.Fs;
    si          = FEData.stft.tp - FEData.stft.wt/2 - FEData.det.wt/2;
    sf          = FEData.stft.tp + FEData.stft.wt/2 + FEData.det.wt/2;
    FEData.stft.tplim  = [si sf];
    
    