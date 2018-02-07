function data = detect(data)
     % Input parameters
        % Spectrogram
        wt      = data.detdata.WindowTime0;
        ovp     = data.detdata.Overlap;
        wf      = data.detdata.WindowFunction;
        Ffilt   = data.detdata.FrequencyFilter;
        % Detection
            Method = data.detdata.MethodKey;
            switch Method
                case 1
                    dfield  = 'Emethod'; 
                    alpha   = data.detdata.Emethod.Threshold;
                    beta    = data.detdata.Emethod.PeakWidth/100;
                    delta   = data.detdata.Emethod.SmoothWindow;
                    theta   = [];
                case 2
                    dfield  = 'VADmethod';
                    alpha   = [];
                    beta    = [];
                    delta   = [];
                    theta = data.detdata.VADmethod.Threshold;
            end
            
    % Routine parameters
    N        = data.audioinfo.FilesAmount;
    Key      = data.audioinfo.Key;    
    SaveMode = data.detdata.SaveMode;
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
    SavePath = cell(N,1);
    
    for i = 1:N
        if Key(i)
            FileTag  = ['D' num2str(data.audioinfo.FileIndex(i),'%06d') ' - ' Name{i}];            
            filepath = fullfile(changepath(data.path.Sonograms,data.path),FileTag);
            SavePath{i} = filepath;
        end
    end  
    
    % Audio signal parameters
    fs      = data.audioinfo.SamplingFrequency;    
    ClipMax = 10;
    w       = nan(N,1);
    ov      = nan(N,1);    
    Smax    = data.audioinfo.Samples;
    
    for i = 1:N
        if Key(i)
            w(i)    = nearest2n (wt,fs(i));
            ov(i)   = round (w(i)*ovp);
        end
    end
    
    %Rutine varibles
    Error                   = cell(N,1);       
    % Spectrogram variables
    
    TimeZero                = zeros(N,1);
    TimeLimits              = zeros(N,2);
    TimeStep                = zeros(N,1);
    FreqZero                = zeros(N,1);
    FreqLimits              = zeros(N,2);
    FreqStep                = zeros(N,1); 
    %ClippingSignal          = cell(N,1);
    
    % Detection variables
    
    DetectionTimeInt        = cell(N,1);
    DetectionTimeLength     = cell(N,1);
    DetectionTimePoint      = cell(N,1);
    
    DetectionFreqInt        = cell(N,1);
    DetectionFreqLength     = cell(N,1);
    DetectionFreqPoint      = cell(N,1);
    
    DetectionPoint          = cell(N,1);
    DetectionSmpPoint       = cell(N,1);
    DetectionPowPoint       = cell(N,1);
    
    DetectionSignal          = cell(N,1);
    NumberOfDetections       = zeros(N,1);        
    DetectionsClipSignal     = cell(N,1);
    
    switch Method
        case 1
            NumberOfHarmonics        = cell(N,1);
        case 2
            DetectionPowRng   = cell(N,1);
            DetectionPowMean  = cell(N,1);
            DetectionLRValue  = cell(N,1); 
            Eccentricity      = cell(N,1);
            Orientation       = cell(N,1);
            Solidity          = cell(N,1);
            MajorAxisLength   = cell(N,1);
            MinorAxisLength   = cell(N,1);
            Area              = cell(N,1);
            FilledArea        = cell(N,1);
            Centroid          = cell(N,1);
            WeightedCentroid  = cell(N,1);
    end
    
    
    
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
            [P,F,T] = spec (y,w(n),ov(n),[],fs(n),wf,'matlab');
        % STFT filtering
            iF = interval(F,Ffilt);
            F0 = max(F(find(iF,1,'first')-1),0);
            f  = F(iF);
            t  = T+T0;
            p  = P(iF,:);

            dT = T(2)-T(1);
            dF = F(2)-F(1);
            
            TimeZero(n)          = T0;
            TimeLimits(n,:)      = [t(1) t(end)];
            TimeStep(n)          = dT;
            FreqZero(n)          = F0;
            FreqLimits(n,:)      = [f(1) f(end)];
            FreqStep(n)          = dF;
            
          % Detection   
            switch Method
                case 1
                    ndelta = round(delta./[dT dF]);
                    [~,Ex,~,~,iDx,~,~,Py,iDyy] = PdetectorXY(p,alpha,beta,ndelta);
                if ~isempty(iDx)
                    Ndet = size(iDx,1);
                    
                    iPt = t(iDx);
                    DetectionTimeInt{n}      = [iPt(:,1) iPt(:,3)];
                    DetectionTimeLength{n}   = iPt(:,3) - iPt(:,1);
                    DetectionTimePoint{n}    = iPt(:,2);
                    

                    iDy = maxiDy (iDyy,Py,Ndet);
                    if Ndet>1
                        iPf = f(iDy);
                    else
                        iPf = f(iDy)';
                    end
                    DetectionFreqInt{n}   = [iPf(:,1) iPf(:,3)];
                    DetectionFreqLength{n}= iPf(:,3) - iPf(:,1);
                    DetectionFreqPoint{n} = iPf(:,2);                   

                    DetectionPoint{n}       = [iDy(:,2) iDx(:,2)];
                    DetectionPowPoint{n}    = p(iDy(:,2), iDx(:,2));
                    DetectionSmpPoint{n}    = round(iPt(:,2)*fs(n));

                    DetectionSignal{n}      = Ex;
                    
                    NumberOfDetections(n)   = Ndet;                
                    if TimeZero(n)==0
                    Clip = isclipped(Y,T,w(n),fs(n),ClipMax);
                    DetectionsClipSignal{n} = nnz(Clip(iDx(:,2)));
                    end
                    
                    Hn = harm(iDyy,f,iPf,Ndet);
                    NumberOfHarmonics{n}    = Hn;
                end                
                case 2
                [vad,ne] = getvadparam (dT,theta);
                [~,LR,~,~,iDx,~,iDy,~,PLR] = VADdetector(p,vad,ne);
                if ~isempty(iDx)
                    Ndet = size(iDx,1);
                    eta  = vad.pr;
                    RefPoint = [iDy' iDx(:,2)];
                    C = callpattern(p,PLR,eta,iDx,RefPoint,Ndet,dF,dT,F0,T0);
                    
                    Tc  = C.BoundingBox(:,1);
                    dTc = C.BoundingBox(:,3);
                    Fc  = C.BoundingBox(:,2);
                    dFc = C.BoundingBox(:,4);
                    
                    DetectionTimeInt{n}      = [Tc Tc+dTc];
                    DetectionTimeLength{n}   = dTc;
                    DetectionTimePoint{n}    = C.MaxIntTime;
                    DetectionSmpPoint{n}     = round(Tc*fs(n));               
                                
                    DetectionFreqInt{n}      = [Fc Fc+dFc];
                    DetectionFreqLength{n}   = dFc;
                    DetectionFreqPoint{n}    = C.MaxIntFreq;               
                
                    DetectionPoint{n}       = RefPoint;
                    DetectionPowPoint{n}    = C.MaxIntensity;
                    DetectionPowRng{n}      = C.MaxIntensity-C.MinIntensity;
                    DetectionPowMean{n}     = C.MeanIntensity;
                    
                    DetectionLRValue{n}     = LR(iDx);
                    DetectionSignal{n}      = LR;
                
                    Eccentricity{n}         = C.Eccentricity;
                    Orientation{n}          = C.Orientation;
                    Solidity{n}             = C.Solidity;
                    MajorAxisLength{n}      = C.MajorAxisLength;
                    MinorAxisLength{n}      = C.MinorAxisLength;
                    Area{n}                 = C.Area;
                    FilledArea{n}           = C.FilledArea;
                    Centroid{n}             = C.Centroid;
                    WeightedCentroid{n}     = C.WeightedCentroid;
                    
                    NumberOfDetections(n)   = Ndet;                
                    if TimeZero(n)==0
                    Clip = isclipped(Y,T,w(n),fs(n),ClipMax);
                    DetectionsClipSignal{n}  = nnz(Clip(iDx(:,2)));
                    savedetvar(p,f,t,C,SavePath{n},SaveMode)
                    end
                end                    
            end
                        
            dtime = tsamples/Smax(n)*toc/poolsize*hxs;
            ftime =  datestr(datenum(itime)+datenum(dtime));
            dispstr = ['Bat call detection. File ' num2str(n) '/' num2str(N) ': ' Name{n} ' - ' num2str(num2str(toc),'%2.1f') 'sec ' ftime ];
            disp(dispstr) 
        catch err
         Key(n) = 0;
         Error{n} = err;
         disp(['Error in file:' num2str(n)]);
         disp(err.message)         
        end       
    end    
end
        if any(~NumberOfDetections)
           Key(~NumberOfDetections) = 0;
           %Error{~NumberOfDetections} = {'Null detections'};
        end
        
        data.detdata.Key                    = Key;
        data.detdata.SavePath               = changepath(SavePath,data.path,computer);           
        data.detdata.Error                  = Error;
        
        data.detdata.WindowSize             = w;        
        data.detdata.WindowTime             = w./fs;
        
        data.detdata.TimeZero   = TimeZero;
        data.detdata.TimeLimits = TimeLimits;
        data.detdata.TimeStep   = TimeStep;
        data.detdata.FreqZero   = FreqZero;
        data.detdata.FreqLimits = FreqLimits;
        data.detdata.FreqStep   = FreqStep;
        
        data.detdata.(dfield).NumberOfDetections     = NumberOfDetections;
        data.detdata.(dfield).DetectionsClipSignal   = DetectionsClipSignal;
        
        data.detdata.(dfield).DetectionTimeInt      = DetectionTimeInt;
        data.detdata.(dfield).DetectionTimeLength   = DetectionTimeLength;
        data.detdata.(dfield).DetectionTimePoint    = DetectionTimePoint;
        
        data.detdata.(dfield).DetectionFreqInt      = DetectionFreqInt;
        data.detdata.(dfield).DetectionFreqLength   = DetectionFreqLength;
        data.detdata.(dfield).DetectionFreqPoint    = DetectionFreqPoint;
        
        data.detdata.(dfield).DetectionPoint        = DetectionPoint;
        data.detdata.(dfield).DetectionPowPoint     = DetectionPowPoint;
        
        data.detdata.(dfield).DetectionSmpPoint     = DetectionSmpPoint;
        data.detdata.(dfield).DetectionSignal       = DetectionSignal;        
        
        switch Method
            case 1
                data.detdata.(dfield).NumberOfHarmonics     = NumberOfHarmonics;
            case 2
                data.detdata.(dfield).DetectionPowRng       = DetectionPowRng;
                data.detdata.(dfield).DetectionPowMean      = DetectionPowMean;
                data.detdata.(dfield).DetectionLRValue      = DetectionLRValue;
                
                data.detdata.(dfield).Eccentricity         = Eccentricity;
                data.detdata.(dfield).Orientation          = Orientation;
                data.detdata.(dfield).Solidity             = Solidity;
                data.detdata.(dfield).MajorAxisLength      = MajorAxisLength;
                data.detdata.(dfield).MinorAxisLength      = MinorAxisLength;
                data.detdata.(dfield).Area                 = Area;
                data.detdata.(dfield).FilledArea           = FilledArea;
                data.detdata.(dfield).Centroid             = Centroid;
                data.detdata.(dfield).WeightedCentroid     = WeightedCentroid;
                data.detdata.(dfield).NumberOfDetections   = NumberOfDetections;                
                data.detdata.(dfield).DetectionsClipSignal = DetectionsClipSignal;
        end

function savedetvar(p,f,t,C,filename,mode)
switch mode
    case 's'
        save(filename,'f','t','p','-v6');
    case 'd'
        save(filename,'C','-v6'); 
    case 's+d'
        save(filename,'f','t','p','C','-v6');
    otherwise            
end                     
function Hn = harm(iDyy,f,iPf,Ndet)
         Hn  = zeros(Ndet,1);
        for i= 1:Ndet
            ipf = f(iDyy{i}(:,2));
            if size(ipf,1)> 1
            ipfmax = iPf(i,2);
            h = false(size(ipf,1),3);
            for j = 1:3
                r = rem(ipf,ipfmax/j)/ (ipfmax/j);
            h(:,j)= 0.9 < r | r < 0.1;                    
            end
            Hn(i) = nnz(sum(h,2));
            end
        end
function iDy = maxiDy (iDyy,Py,Ndet)
        iDy = zeros(Ndet,3);
        for i = 1:Ndet
            idy = iDyy{i};
            if size(idy,1)> 1
                ePy = Py{i};
                [~,j] = max(ePy);                        
                iDy(i,:) = idy(j,:);
            else
                iDy(i,:) = idy;                        
            end
        end
function [vad,ne] = getvadparam (dT,theta)
            vad.dt  = dT;             % true frame increment time
            vad.pr  = theta;            % Speech probability threshold
            vad.gx  = db2pow(30);     % maximum posterior SNR = 30dB
            vad.gz  = db2pow(-40);    % minimum posterior SNR = -40dB

            vad.ne  = 0;              % noise estimation: 0=min statistics, 1=MMSE [0]

            vad.ge  = 1;                % xi estimation: 0= Itakura Saito (ISD), 1 = Decision Directed (DD) [1]
            vad.ta  = -dT/log(0.98);    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [2]
            vad.xn  = 0;                % minimum prior SNR = -Inf dB

            vad.hmm = 1;                % HMM-Based Hang Over [0]
            vad.ts  = 0.04;             % mean talkspurt length (100 ms)
            vad.tn  = 0.005;             % mean silence length (50 ms)


            % Estimate noise spectrum using minimum statistics
            ne.taca   =-dT/log(0.95);       % smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11)
            ne.tamax  =-dT/log(0.96);      % max smoothing time constant in (3) = -tinc/log(0.96)
            ne.taminh =-dT/log(0.3);       % min smoothing time constant (upper limit) in (3) = -tinc/log(0.3)
            ne.tpfall =0.064/10;           % time constant for P to fall (12)
            ne.tbmax  =-dT/log(0.8);             % max smoothing time constant in (20) = -tinc/log(0.8)
            ne.qeqmin =2;                  % minimum value of Qeq (23)
            ne.qeqmax =14;                 % max value of Qeq per frame
            ne.av     =2.12;               % fudge factor for bc calculation (23 + 13 lines)
            ne.td     =0.01;              % time to take minimum over
            ne.nu     =16;                  % number of subwindows
            ne.qith   =[0.03 0.05 0.06 Inf];% noise slope thresholds in dB/s
            ne.nsmdb  =[47 31.4 15.7 4.1];    
function  C = callpattern(p,PLR,eta,iPx,RefPoint,N,dF,dT,F0,T0)
 
    [I,J] = size(p);
    Deta = PLR >=eta;
    
    iRefPoint = sub2ind([I J],RefPoint(:,1),RefPoint(:,2));
    
    C = struct();
    C.Connectivity  = 8;
    C.ImageSize     = [I,J];
    C.NumObjects    = N;
    C.PixelIdxList  = cell(1,N);
   
    for i = 1:N
        Dc = false(I,J);
        j = iPx(i,1):iPx(i,3);
        Dc(:,j) = 1;
        BW = Deta & Dc;
        c = bwconncomp(BW,8);
        for k = 1:c.NumObjects
            if any(c.PixelIdxList{k} == iRefPoint(i))
                    C.PixelIdxList{i}= c.PixelIdxList{k};
                    C.RefPoint(i,:)  = RefPoint(i,:);
            end
        end
    end
    
    Cprops = regionprops(C,p,'Area','FilledArea','Centroid','BoundingBox',...
                            'Eccentricity','FilledArea','Orientation','Solidity',...
                            'MajorAxisLength','MinorAxisLength','WeightedCentroid',...
                            'MaxIntensity','MeanIntensity','MinIntensity',...
                            'PixelList','PixelValues');
    B=C;                    
    fnames = fieldnames(Cprops);                
    for i = 1: length(fnames)
        switch fnames{i}
            case {'PixelList','PixelValues'}
                B.(fnames{i}) = cat(2,{Cprops.(fnames{i})});
            otherwise
                B.(fnames{i}) = cat(1,Cprops.(fnames{i}));
        end
    end
    
    C.Eccentricity      = B.Eccentricity;
    C.Orientation       = B.Orientation;
    C.Solidity          = B.Solidity;
    C.MajorAxisLength   = B.MajorAxisLength;
    C.MinorAxisLength   = B.MinorAxisLength;
    C.Area              = B.Area*dF*dT;
    C.FilledArea        = B.FilledArea*dF*dT;
    C.Centroid          = [B.Centroid(:,1)*dT+T0          B.Centroid(:,2)*dF+F0];
    C.WeightedCentroid  = [B.WeightedCentroid(:,1)*dT+T0  B.WeightedCentroid(:,2)*dF+F0];
    C.BoundingBox       = B.BoundingBox.*repmat([dT dF dT dF],C.NumObjects,1) + repmat([T0 F0 0 0],C.NumObjects,1);
    C.MaxIntensity      = pow2db(B.MaxIntensity);
    C.MeanIntensity     = pow2db(B.MeanIntensity);
    C.MinIntensity      = pow2db(B.MinIntensity);
    C.PixelList         = B.PixelList;
    C.PixelValues       = B.PixelValues;
    
    for i =1:C.NumObjects
        ip = find(B.PixelValues{i} == B.MaxIntensity(i),1,'first');
    C.MaxIntFreq(i,1)   = B.PixelList{i}(ip,2)*dF+F0;
    C.MaxIntTime(i,1)   = B.PixelList{i}(ip,1)*dT+T0;
    end
    
    