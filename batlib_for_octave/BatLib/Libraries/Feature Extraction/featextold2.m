function data = featextx(data)

    det     = detparameters(data);    
    psd     = psdparameters(data);
    stft    = stftparameters(data);   
    
    data.featext.PSD.Pxx                = cell(data.audioinfo.FilesAmount,1);
    data.featext.PSD.Frequency          = cell(data.audioinfo.FilesAmount,1);
    data.featext.PSD.PeakFrequency      = cell(data.audioinfo.FilesAmount,1);
    data.featext.PSD.PeakValue          = cell(data.audioinfo.FilesAmount,1);
    data.featext.PSD.PeakFrequencyMean  =zeros(data.audioinfo.FilesAmount,1);
    data.featext.PSD.PeakFrequencyStd   =zeros(data.audioinfo.FilesAmount,1);
    
    data.featext.STFT.EFTcurve          = cell(data.audioinfo.FilesAmount,1);
    data.featext.STFT.InitialTime       = cell(data.audioinfo.FilesAmount,1);
for n = find(data.audioinfo.Key)'
    disp(['Feature extraction. File ' num2str(n) ': ' data.audioinfo.Name{n}]) 
    %Audio parameters
    path = data.audioinfo.Path{n};
    fs   = data.audioinfo.SamplingFrequency(n);
    %Detection parameters
    det = detparameters(data,det,fs,n);        
    %PSD parameters
    psd  = psdparameters(data,psd,fs);
    %STFT parameters
    stft = stftparameters(data,stft,fs);
    
    %Selection of detection
    i       = data.detdata.Key{n}; 
    M       = data.detdata.SelectionAmount(n);
    % Samples limits
    [det stft] = getsamplelimits(det,stft,fs,i);
    
    PSD     = zeros(M,psd.nfft/2+1);
    Fp      = zeros(M,1);
    maxPSD  = zeros(M,1);
    
    EFTcurve   = cell(M,3);
    EFTt0      = zeros(M,1);
    
    EFTfeat    = 0;
    for m = 1:M
    disp(['Feature extraction. File ' num2str(n) ': Call ' num2str(m) ' of ' num2str(M)]) 
        %Audio reading
            samples = det.S0(m,:);
            audio   = wavread(path,samples);
        %PSD features
            samples = det.s(m,1):det.s(m,2);
            Xpsd = psdfeatures (audio(samples),psd,fs);
            PSD(m,:)    = Xpsd{1};
            Fpsd        = Xpsd{2};
            maxPSD(m)   = Xpsd{3};
            Fp(m)       = Xpsd{4};
        %STFT features
            %Spectrogram frame
                [P,F,T,i0,j0,t0,Pi,Fi,Ti] = stftframe(audio,stft,fs,m);
            %Get signal distribution
                S = getsignaldist(P,i0,j0,stft);
            %EFT Contour
                [e,f,t] = EFTcontour (P,F,T,S);
                EFTcurve{m,1} = e;
                EFTcurve{m,2} = f;
                EFTcurve{m,3} = t;
                EFTt0(m)      = t0 + det.S0(m,1)/fs;
                
                Y = EFTfeatures ([],e,f,t);
                %plotcall(P,F,T,S,Pi,Fi,Ti,i0,j0,e,f,t,stft,m)
                if any(isnan(e))|| any(isnan(f))|| any(isnan(t))
                    plotcall(P,F,T,S,Pi,Fi,Ti,i0,j0,e,f,t,stft,m)
                    disp(find(isnan(e)));
                end

     end

    data.featext.PSD.Pxx{n}                 = PSD;
    data.featext.PSD.Frequency{n}           = Fpsd;
    data.featext.PSD.PeakFrequency{n}       = Fp;
    data.featext.PSD.PeakValue{n}           = maxPSD;
    data.featext.PSD.PeakFrequencyMean(n)   = mean(Fp);
    data.featext.PSD.PeakFrequencyStd(n)    = std(Fp);
    
    data.featext.STFT.EFTcurve{n}           = EFTcurve;
    data.featext.STFT.InitialTime{n}        = EFTt0;
end