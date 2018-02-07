function X = psdfeatures (audio,psdp,fs,m)
if nargin == 0
    X = {'Pxx','Frequency','PeakValue','PeakFrequency'};
else
        %PSD
        [P,F] = pwelch(audio,psdp.w,psdp.ov,psdp.nfft,fs);
        iFfilt = interval(F,psdp.flim(m,:));
        %Peak Frequency
        [maxP,iFp0] = max(P(iFfilt)); 
        iFp = iFp0 + find(iFfilt,1,'first')-1;
        maxF = F(iFp);
    X = {P,F,maxP,maxF};
end
