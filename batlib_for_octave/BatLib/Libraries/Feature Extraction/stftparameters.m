function stft = stftparameters(data,i)
    
    stft.wt       = data.featext.STFT.WindowTime0;
    stft.ovp      = data.featext.STFT.Overlap;
    stft.wf       = data.featext.STFT.WindowType;
    %eval(['stft.whandle = @' stft.wtype ';'])
    stft.df       = data.featext.STFT.deltaF;
    stft.dt       = data.featext.STFT.deltaT;
    stft.aframe   = stft.df*stft.dt;
    stft.amax     = data.featext.EFT.MaximumArea/stft.aframe;
    stft.amin     = data.featext.EFT.MinimumArea/stft.aframe;    
    stft.thetha   = data.featext.EFT.Threshold;
    stft.delta    = data.featext.EFT.AdjustmentStep;
    stft.mode     = data.featext.EFT.FollowFreqPeak;
    stft.conn     = data.featext.EFT.Connectivity;
    stft.areatype = data.featext.EFT.AreaType;
    
    fs = data.audioinfo.SamplingFrequency(i);
    wn = nearest2n (stft.wt,fs);
    stft.fs     = fs; 
    stft.wn     = wn;
    stft.ov     = round (wn*stft.ovp);        
    stft.dS     = round(stft.dt*fs);
    stft.fstep  = fs/wn;
    stft.tstep  = (wn-stft.ov)/fs;   
end