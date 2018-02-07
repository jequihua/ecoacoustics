function data = detect(data)

    N        = data.audioinfo.FilesAmount;
    Key      = data.audioinfo.Key;
    Name     = data.audioinfo.Name;
    SavePath = data.detdata.SavePath;
    System   = data.path.Sytem;
    WinHeader= data.path.WinHeader;
    LinHeader= data.path.LinHeader;    
    Error    = data.detdata.Error;
    
    %wt                 = data.detdata.WindowTime;
    SamplingFrequency  = data.audioinfo.SamplingFrequency;
    %SamplesLength      = data.audioinfo.Samples;
    TimeStep           = data.detdata.TimeStep;
    FreqStep           = data.detdata.FreqStep;
    TimeZero           = data.detdata.TimeZero;
    FreqZero           = data.detdata.FreqZero;
    ClippingSignal     = data.detdata.ClippingSignal;
    
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
    %DetectionSignalPeakIndex = cell(N,1);
    NumberOfDetections       = zeros(N,1);
    NumberOfHarmonics        = cell(N,1);    
    DetectionsClipSignal     = cell(N,1);
    
    try
        Pool                    = gcp('nocreate');
        if isempty(Pool)
        poolsize = 1;
        else
        poolsize = Pool.NumWorkers;
        end
    catch
        poolsize=12;
    end
    
switch data.detdata.MethodKey
    case 1
    alpha   = data.detdata.Emethod.Threshold;
    beta    = data.detdata.Emethod.PeakWidth/100;
    delta   = data.detdata.Emethod.SmoothWindow;
    TFstep  = [TimeStep FreqStep];
    ndelta  = round(repmat(delta,size(TFstep,1),1)./TFstep);
    
    
    itime = now;
    hxs = 1/(60*60*24);
    for n = 1:N
        if Key(n)
            try
            tic
            
            filepath = changepath(SavePath{n},System,WinHeader,LinHeader); 
            [p,f,t]  =  readsonogram(filepath); 
            
            [~,Ex,~,Px,iDx,Ey,~,Py,iDyy] = PdetectorXY(p,alpha,beta,ndelta(n,:));
            
                if ~isempty(iDx)
                
                Fs   = SamplingFrequency(n);
                Ndet = size(iDx,1);                
                Clip = ClippingSignal{n};
                
                iPt = t(iDx);
                DetectionTimeInt{n}      = [iPt(:,1) iPt(:,3)];
                DetectionTimeLength{n}   = iPt(:,3) - iPt(:,1);
                DetectionTimePoint{n}    = iPt(:,2);
                DetectionSmpPoint{n}     = round(iPt(:,2)*Fs);
                
                iDy = maxiDy (iDyy,Py,Ndet);
                if Ndet>1
                    iPf = f(iDy);
                else
                    iPf = f(iDy)';
                end
                DetectionFreqInt{n}   = [iPf(:,1) iPf(:,3)];
                DetectionFreqLength{n}= iPf(:,3) - iPf(:,1);
                DetectionFreqPoint{n} = iPf(:,2);
                
                Hn = harm(iDyy,f,iPf,Ndet);
                NumberOfHarmonics{n}    = Hn;
                
                DetectionPoint{n}       = [iDy(:,2) iDx(:,2)];
                DetectionPowPoint{n}    = p([iDy(:,2) iDx(:,2)]);
                
                NumberOfDetections(n)   = Ndet;                
                if TimeZero(n)==0                
                DetectionsClipSignal{n}  = nnz(Clip(iDx(:,2)));
                end
                
                % Disp Info
                ftime =  datestr(datenum(itime)+datenum(N*toc*hxs/poolsize));
                dispstr = ['Detection. File ' num2str(n) '/' num2str(N) ': ' Name{n} ' - ' num2str(num2str(toc),'%2.1f') 'sec ' ftime ' N: ' num2str(Ndet,'%1.0d')  ];
                disp(dispstr)
                end
                
                catch err
                Key(n) = 0;
                Error{n} = err.message;
                disp(err.message)
            end
        end
    end

    data.detdata.Emethod.DetectionTimeInt        = DetectionTimeInt;
    data.detdata.Emethod.DetectionTimeLength     = DetectionTimeLength;
    data.detdata.Emethod.DetectionTimePoint      = DetectionTimePoint;
    
    data.detdata.Emethod.DetectionFreqInt        = DetectionFreqInt;
    data.detdata.Emethod.DetectionFreqLength     = DetectionFreqLength;
    data.detdata.Emethod.DetectionFreqPoint      = DetectionFreqPoint;
    
    data.detdata.Emethod.DetectionPoint          = DetectionPoint;
    data.detdata.Emethod.DetectionSmpPoint       = DetectionSmpPoint;
    data.detdata.Emethod.DetectionPowPoint       = DetectionPowPoint;
    
    %DetectionSignal          = cell(N,1);
    %DetectionSignalPeakIndex = cell(N,1);
    data.detdata.Emethod.NumberOfDetections       = NumberOfDetections;
    data.detdata.Emethod.NumberOfHarmonics        = NumberOfHarmonics;    
    data.detdata.Emethod.DetectionsClipSignal     = DetectionsClipSignal;
    
    data.detdata.Key                    = Key;
    data.detdata.Error                  = Error;
    
    
    case 2
        
    theta = data.detdata.VADmethod.Threshold;
    itime = now;
    hxs = 1/(60*60*24);
    
    DetectionPowRng   = cell(N,1);
    DetectionPowMean  = cell(N,1); 
    Eccentricity      = cell(N,1);
    Orientation       = cell(N,1);
    Solidity          = cell(N,1);
    MajorAxisLength   = cell(N,1);
    MinorAxisLength   = cell(N,1);
    Area              = cell(N,1);
    FilledArea        = cell(N,1);
    Centroid          = cell(N,1);
    WeightedCentroid  = cell(N,1);
    
    SavePathC         = cell(N,1);
    for i = 1:N
        if Key(i)
            FileTag  = ['C' num2str(data.audioinfo.FileIndex(i),'%06d') ' - ' Name{i}];            
            filepath = fullfile( changepath(data.path.Sonograms,System,WinHeader,LinHeader),FileTag);
            SavePathC{i} = filepath;
        end
    end
    
    parfor n = 1:N
        if Key(n)
            try
            tic
            
            filepath = changepath(SavePath{n},System,WinHeader,LinHeader); 
            [p,~,~]  = readsonogram(filepath);
            [vad,ne] = getvadparam (TimeStep(n),theta);
            [~,LR,~,~,iDx,~,iDy,~,PLR] = VADdetector(p,vad,ne);
            if ~isempty(iDx)
                
                Fs   = SamplingFrequency(n);
                Ndet = size(iDx,1);                
                Clip = ClippingSignal{n};
                eta  = vad.pr;
                RefPoint = [iDy' iDx(:,2)];
                C = callpattern(p,PLR,eta,iDx,RefPoint,Ndet,FreqStep(n),TimeStep(n),FreqZero(n),TimeZero(n));
                savecallpattern(C,SavePathC{n})
                
                T  = C.BoundingBox(:,1);
                dT = C.BoundingBox(:,3);
                F  = C.BoundingBox(:,2);
                dF = C.BoundingBox(:,4);
                
                DetectionTimeInt{n}      = [T T+dT];
                DetectionTimeLength{n}   = dT;
                DetectionTimePoint{n}    = C.MaxIntTime;
                DetectionSmpPoint{n}     = round(T*Fs);               
                                
                DetectionFreqInt{n}      = [F F+dF];
                DetectionFreqLength{n}   = dF;
                DetectionFreqPoint{n}    = C.MaxIntFreq;               
                
                DetectionPoint{n}       = RefPoint;
                DetectionPowPoint{n}    = C.MaxIntensity;
                DetectionPowRng{n}      = C.MaxIntensity-C.MinIntensity;
                DetectionPowMean{n}     = C.MeanIntensity;
                DetectionSignal{n}      = LR;
                
                Eccentricity{n}      = C.Eccentricity;
                Orientation{n}       = C.Orientation;
                Solidity{n}          = C.Solidity;
                MajorAxisLength{n}   = C.MajorAxisLength;
                MinorAxisLength{n}   = C.MinorAxisLength;
                Area{n}              = C.Area;
                FilledArea{n}        = C.FilledArea;
                Centroid{n}          = C.Centroid;
                WeightedCentroid{n}  = C.WeightedCentroid;
                NumberOfDetections(n)   = Ndet;                
                if TimeZero(n)==0                
                DetectionsClipSignal{n}  = nnz(Clip(iDx(:,2)));
                end
            % Disp Info
            ftime =  datestr(datenum(itime)+datenum(N*toc*hxs/poolsize));
            dispstr = ['Detection. File ' num2str(n) '/' num2str(N) ': ' Name{n} ' - ' num2str(num2str(toc),'%2.1f') 'sec ' ftime ' N: ' num2str(Ndet,'%1.0d')  ];
            disp(dispstr)
            end
            
            catch err
                Key(n) = 0;
                Error{n} = err.message;
                disp(err.message)
            end
        end
    end
    
                data.detdata.VADmethod.DetectionTimeInt      = DetectionTimeInt;
                data.detdata.VADmethod.DetectionTimeLength   = DetectionTimeLength;
                data.detdata.VADmethod.DetectionTimePoint    = DetectionTimePoint;
                data.detdata.VADmethod.DetectionSmpPoint     = DetectionSmpPoint;
                data.detdata.VADmethod.DetectionFreqInt      = DetectionFreqInt;
                data.detdata.VADmethod.DetectionFreqLength   = DetectionFreqLength;
                data.detdata.VADmethod.DetectionFreqPoint    = DetectionFreqPoint;
                data.detdata.VADmethod.DetectionPoint        = DetectionPoint;
                data.detdata.VADmethod.DetectionPowPoint     = DetectionPowPoint;
                data.detdata.VADmethod.DetectionPowRng       = DetectionPowRng;
                data.detdata.VADmethod.DetectionPowMean      = DetectionPowMean;
                data.detdata.VADmethod.DetectionSignal       = DetectionSignal;
                
                data.detdata.VADmethod.Eccentricity         = Eccentricity;
                data.detdata.VADmethod.Orientation          = Orientation;
                data.detdata.VADmethod.Solidity             = Solidity;
                data.detdata.VADmethod.MajorAxisLength      = MajorAxisLength;
                data.detdata.VADmethod.MinorAxisLength      = MinorAxisLength;
                data.detdata.VADmethod.Area                 = Area;
                data.detdata.VADmethod.FilledArea           = FilledArea;
                data.detdata.VADmethod.Centroid             = Centroid;
                data.detdata.VADmethod.WeightedCentroid     = WeightedCentroid;
                data.detdata.VADmethod.NumberOfDetections   = NumberOfDetections;                
                data.detdata.VADmethod.DetectionsClipSignal = DetectionsClipSignal;
                
                data.detdata.Key                    = Key;
                data.detdata.Error                  = Error;
        
end
function [p,f,t] =  readsonogram(filepath)
            sonogram  = load(filepath);
            p  = sonogram.p;    
            t  = sonogram.t;
            f  = sonogram.f;
            
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
    
      function savecallpattern(C,filename)
          save(filename,'C','-v6');
            
