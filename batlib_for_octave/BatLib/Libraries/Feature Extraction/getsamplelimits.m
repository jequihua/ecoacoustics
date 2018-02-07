function [det,stft] = getsamplelimits(det,stft,i)

    fs = det.fs;

    Sp = det.Sp(i);
    Si = det.S(i,1);
    Sf = det.S(i,2);

    S0i = round(Sp - max(Sp-Si,stft.dS/2) - det.dSp/2);
    S0i = max(S0i,1);
    S0f = round(Sp + max(Sf-Sp,stft.dS/2) + det.dSp/2);
    S0f = min(S0f,det.Smax);
    
    det.S0 = [S0i S0f];
    det.s  = [Si Sf] - [S0i S0i] +1;
    det.T0 = det.S0(:,1)/fs;

    stft.sp     = Sp - S0i +1;
    stft.tp     = stft.sp/fs;
    
    si          = stft.tp - stft.wt/2 - det.wt/2;
    sf          = stft.tp + stft.wt/2 + det.wt/2;
    stft.tplim  = [si sf];
    
    fi          = det.dFp(i,1) - 2/stft.wt - 2/det.wt;
    ff          = det.dFp(i,2) + 2/stft.wt + 2/det.wt;
    stft.fplim  = [fi ff];


