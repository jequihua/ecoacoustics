function data = CricketFilter(data,obj)
     % Input parameters
     
        % obj
         
         CrClass = obj.RFclass; 
         D       = obj.D;
         CFilt = obj.CFilt;
         Tw    = obj.UndSampTime;
         Fw    = obj.UndSampFreq;
         nTw   = obj.UndSampAmount;
         nMus  = obj.PMUSICdim;
         
         nPFpks   = obj.PFNumPeaks;
         hPFpmin  = obj.PFPeakHeigthMin;
         wPFp     = obj.PFPeakWeightLevel;
         hPTpmin   = obj.TMinPeakHeigth;
         
         
        % Spectrogram
        wt      = data.detdata.WindowTime;
        wn      = data.detdata.WindowSize;
        ovp     = data.detdata.Overlap;
        wf      = data.detdata.WindowFunction;
        Ffilt   = data.detdata.FrequencyFilter;
        % Detection
            Method = data.detdata.MethodKey;
            switch Method
                case 1
                    dfield  = 'Emethod'; 
                case 2
                    dfield  = 'VADmethod';                    
            end
            
    % Routine parameters
    N        = data.audioinfo.FilesAmount;
    Key      = data.audioinfo.Key;    
        try
            Pool                    = gcp('nocreate');
            if isempty(Pool)
            poolsize = 1;
            else
            poolsize = Pool.NumWorkers
            end
        catch
            poolsize=4;
        end
    tsamples = nansum(data.audioinfo.Samples(Key));    
    itime = now;
    hxs = 1/(60*60*24);
    
    % Audio file info 
    
    System   = data.path.Sytem;
    WinHeader= data.path.WinHeader;
    LinHeader= data.path.LinHeader;
    
    Name     = data.audioinfo.Name;
    Path     = changepath(data.audioinfo.Path,data.path);
    
    % Audio signal parameters
    fs      = data.audioinfo.SamplingFrequency;    
    Smax    = data.audioinfo.Samples;    
    for i = 1:N
        if Key(i)
           ov(i)   = round (wn(i)*ovp);
        end
    end
    
    %Rutine varibles
    Error                   = cell(N,1);       
    % Spectrogram variables
    
    AudioFeatures           = zeros(N,D);
    CricketKey              = false(N,1);
    CricketProb             = zeros(N,1);
    
    % Detection variables
    
    NumberOfDetections              = data.detdata.(dfield).NumberOfDetections;        
    DetectionTimeInt                = data.detdata.(dfield).DetectionTimeInt;
    DetectionFreqInt                = data.detdata.(dfield).DetectionFreqInt;
    
    
    CricketDetectionTimeInt        = cell(N,1);
    CricketDetectionTimeLength     = cell(N,1);
    CricketDetectionTimePoint      = cell(N,1);
    
    CricketDetectionFreqInt        = cell(N,1);
    CricketDetectionFreqLength     = cell(N,1);
    CricketDetectionFreqPoint      = cell(N,1);
    
    CricketSelectionKey            = cell(N,1);
    
parfor n = 1:N
    if Key(n)
        try
            tic              
        %Audio reading
            Y     = audioread(Path{n});
            Y     = Y(:,1); 
            S0    = find(Y,1,'first');
            T0    = (S0-1)/fs(n);    
            y = Y(S0:end)/max(Y)*0.9;        
        % SFTF calculation
            [P,F,T] = spec (y,wn(n),ov(n),[],fs(n),wf,'tftb');
        % STFT filtering
            iF = interval(F,Ffilt);
            F0 = max(F(find(iF,1,'first')-1),0);
            f  = F(iF);
            t  = T+T0;
            p  = P(iF,:);

            dT = T(2)-T(1);
            dF = F(2)-F(1);
            
            PdB   = pow2db(p);
            Pstd  = repmat(std(PdB,[],2),1,size(PdB,2));
            P0    = repmat(median(PdB,2),1,size(PdB,2));
            P     = (PdB-P0)./Pstd; 
            PF    = mean(P,2);
            
            iFc = interval(f,CFilt);
            PT  = mean(p(iFc,:),1);
            
            
            % Feature extraction
                Fpf  = zeros(nPFpks,1);
                Epf  = zeros(nPFpks,1);
                Wpf  = zeros(nPFpks,1);
                Hpf  = zeros(nPFpks,1);
                
                Fpt  = zeros(nPFpks,1);
                Ept  = zeros(nPFpks,1);
                Wpt  = zeros(nPFpks,1);
                Hpt  = zeros(nPFpks,1);
            
                % PF space
                
                [pks,ipks,hpks,wpks,intpks] = peakdet2(PF,hPFpmin,nPFpks,'descend','height',wPFp); 
                npks = length(ipks);
                Fpf(1:npks,1)  = f(ipks);
                Epf(1:npks,1)  = pks;
                Wpf(1:npks,1)  = wpks*dF;
                Hpf(1:npks,1)  = hpks;
                Ipf            = f(intpks);
                
                % PT space
                pt      = smooth(PT,nTw);
                pt      = decimate(pt,nTw);
                pt      = (pt-mean(pt))/std(pt);
                [Pxx,Fxx] = pmusic(pt,nMus,[],Fw); 
                Pxx = pow2db(Pxx);
                dFxx = Fxx(2)-Fxx(1);
                
                [pks,ipks,hpks,wpks] = peakdet2(Pxx,hPTpmin,nMus,'descend');
                npks = length(ipks);
                Fpt(1:npks,1)  = Fxx(ipks);
                Ept(1:npks,1)  = pks;
                Wpt(1:npks,1)  = wpks*dFxx;
                Hpt(1:npks,1)  = hpks;
            
                 x      = [Fpf' Epf' Wpf' Hpf' Fpt' Ept' Wpt' Hpt'];
                 [y,prob]  =  predict(CrClass,x);
                 y      = logical(str2double(y));
                
                    AudioFeatures (n,:)     = x;
                    CricketKey(n)           = y;
                    CricketProb(n)          = max(prob);
                    Ndet  = NumberOfDetections(n);
                    iC    = false(Ndet,1);
                 if y
                     disp(['Cricket detected File ' num2str(n) '/' num2str(N) ' Freq:' num2str(Fpf(1)*1e-3,'%2.0f kHz') ' Pulse Int:' num2str(1/Fpt(1),'%2.2f sec') ]);
                     if interval(Fpf(1),Ipf(1,:)) 
                        iFc    = interval(f,Ipf(1,:));
                     else
                         iFc  = interval(f,CFilt);
                     end
                     pt  = mean(p(iFc,:),1);
                     pt  = pt./std(pt); 
                     pt = smooth(pt,nTw);
                    [pks,ipks,~,wpks,intpks] = peakdet2(pt,0.5);
                    npks  = length(ipks);
                    dfpks = [min(f(iFc)) max(f(iFc))];
                    dtpks =  t(intpks);
                    
                    dT    = DetectionTimeInt{n};
                    dF    = DetectionFreqInt{n};                  
                    
                    for i = 1:Ndet
                        for j = 1: npks
                            if dtpks(j,1) < dT(i,2) &  dT(i,1) < dtpks(j,2) 
                                if dfpks(1) < dF(i,2) &  dF(i,1) < dfpks(2)
                                    iC(i) = 1;
                                end
                            end
                        end
                    end
                    disp(['Calls deselected: ' num2str(nnz(iC))]);
                    CricketSelectionKey{n}          = iC;
                    
                    CricketDetectionFreqInt{n}      = dfpks;
                    CricketDetectionFreqLength{n}   = dfpks(2) - dfpks(1);
                    CricketDetectionFreqPoint{n}    = Fpf(1);
                    
                    CricketDetectionTimeInt{n}      = dtpks;
                    CricketDetectionTimeLength{n}   = t(intpks);
                    CricketDetectionTimePoint{n}    = t(ipks);
                 end      

                        
            dtime = tsamples/Smax(n)*toc/poolsize*hxs;
            ftime =  datestr(datenum(itime)+datenum(dtime));
            dispstr = ['Cricket detection. File ' num2str(n) '/' num2str(N) ': ' Name{n} ' - ' num2str(num2str(toc),'%2.1f') 'sec ' ftime ];
            disp(dispstr) 
        catch err
             Key(n) = 0;
             Error{n} = err;
             disp(['Error in file:' num2str(n)]);
             disp(err.getReport)         
         end
    end
end    

        data.detdata.CricketFilter.Key                   = CricketKey;
        data.detdata.CricketFilter.Error                 = Error;
        data.detdata.CricketFilter.AudioFeatures         = AudioFeatures;
        data.detdata.CricketFilter.DetectionProb         = CricketProb;
        data.detdata.CricketFilter.SelectionKey          = CricketSelectionKey;
        
        data.detdata.CricketFilter.DetectionTimeInt      = CricketDetectionTimeInt;
        data.detdata.CricketFilter.DetectionTimeLength   = CricketDetectionTimeLength;
        data.detdata.CricketFilter.DetectionTimePoint    = CricketDetectionTimePoint;
        
        data.detdata.CricketFilter.DetectionFreqInt      = CricketDetectionFreqInt;
        data.detdata.CricketFilter.DetectionFreqLength   = CricketDetectionFreqLength;
        data.detdata.CricketFilter.DetectionFreqPoint    = CricketDetectionFreqPoint;