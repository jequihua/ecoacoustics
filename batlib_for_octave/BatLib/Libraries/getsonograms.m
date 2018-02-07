function data = getsonograms(data)
    wt      = data.detdata.WindowTime;
    ovp     = data.detdata.Overlap;
    Ffilt   = data.detdata.FrequencyFilter;
    
    N    = data.audioinfo.FilesAmount;
    Key = data.audioinfo.Key;
    
    Name = data.audioinfo.Name;
    System                  = data.path.Sytem;
    WinHeader               = data.path.WinHeader;
    LinHeader               = data.path.LinHeader;
    Path = changepath(data.audioinfo.Path,System,WinHeader,LinHeader);
    SamplingFrequency   = data.audioinfo.SamplingFrequency;
    
    ClipMax = 10;
    

    TimeZero                = zeros(N,1);
    TimeLimits              = zeros(N,2);
    TimeStep                = zeros(N,1);
    FreqZero                = zeros(N,1);
    FreqLimits              = zeros(N,2);
    FreqStep                = zeros(N,1); 
    ClippingSignal          = cell(N,1);    
    
    SavePath                = cell(N,1);
    
    Error                   = cell(N,1);    
    Pool                    = gcp('nocreate');
    if isempty(Pool)
    poolsize = 1;
    else
    poolsize = Pool.NumWorkers;
    end
    for i = 1:N
        if Key(i)
            FileTag  = ['F' num2str(data.audioinfo.FileIndex(i),'%06d') ' - ' Name{i}];            
            filepath = fullfile( changepath(data.path.Sonograms,System,WinHeader,LinHeader),FileTag);
            SavePath{i} = filepath;
        end
    end
    
    itime = now;
    hxs = 1/(60*60*24);  
    
parfor n = 1:N
    if Key(n)
        try
            tic
        
        %Audio parameters
            path = Path{n};
            fs   = SamplingFrequency(n);    
        %STFT parameters
            w = nearest2n (wt,fs);
            ov = round (w*ovp);        
        %Audio reading
            Y     = audioread(path);
            Y     = Y(:,1); 
            S0    = find(Y,1,'first');
            T0    = (S0-1)/fs;    
            y = Y(S0:end)/max(Y)*0.9;        
        % SFTF calculation
            [~,F,T,P] = spectrogram (y,w,ov,[],fs);
        % Filtering
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

         % Clipping Test
            ClippingSignal{n} = isclipped(Y,T,w,fs,ClipMax);

         % SaveData            
            savespec(t,f,p,SavePath{n});
         % Disp Info
            
            ftime =  datestr(datenum(itime)+datenum(N*toc*hxs/poolsize));
            dispstr = ['Calculating Sonogram. File ' num2str(n) '/' num2str(N) ': ' Name{n} ' - ' num2str(num2str(toc),'%2.1f') 'sec ' ftime ]
            disp(dispstr) 
        catch err
         Key(n) = 0;
         Error{n} = err.message;
         disp(err.message)
        end       
    end
end

    data.detdata.TimeZero   = TimeZero;
    data.detdata.TimeLimits = TimeLimits;
    data.detdata.TimeStep   = TimeStep;
    data.detdata.FreqZero   = FreqZero;
    data.detdata.FreqLimits = FreqLimits;
    data.detdata.FreqStep   = FreqStep;    
    data.detdata.ClippingSignal = ClippingSignal;
    
    data.audioinfo.Key                  = Key;
    data.detdata.SavePath               = SavePath;
    %data.detdata.DetectionTimeLength    = DetectionTimeLength;    
    data.detdata.Error                  = Error;
    
end
function savespec(t,f,p,filename)
        save(filename,'f','t','p','-v6');
    end