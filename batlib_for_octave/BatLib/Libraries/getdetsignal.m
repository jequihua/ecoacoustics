function data = getdetsignal(data)
    wt      = data.detdata.WindowTime;
    ovp     = data.detdata.Overlap;
    Ffilt   = data.detdata.FrequencyFilter;
    
    N    = data.audioinfo.FilesAmount;
    Key = data.audioinfo.Key;
    
    Name = data.audioinfo.Name;
    Path = data.audioinfo.Path;
    SamplingFrequency   = data.audioinfo.SamplingFrequency;
    
    ClipMax = 10;
    

    EnergySignal            = cell(N,1);
    EnergySignaldB          = cell(N,1);
    TimeSignal              = cell(N,1);
    TimeStep                = zeros(N,1);
    ClippingKey             = cell(N,1);    
    
    SavePath                = cell(N,1);    
    Error                   = cell(N,1);
    FileName                = cell(N,1);
    for i = 1:N
        if Key(i)
            FileTag  = ['F' num2str(data.audioinfo.FileIndex(i),'%06d') ' - ' Name{i}];  
            FileName{i} = fullfile(data.path.Sonograms,FileTag);
        end
    end
    
   
    
parfor n = 1:N
    if Key(n)
        try
    disp(['Detection signal. File ' num2str(n) ': ' Name{n}]) 
    %Audio parameters
        path = Path{n};
        fs   = SamplingFrequency(n);    
    %STFT parameters
        w = nearest2n (wt,fs);
        ov = round (w*ovp);        
    %Audio reading
        Y     = audioread(path);
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
        
     % Clipping Test
        ClippingSignal{n} = isclipped(Y,T,w,fs,ClipMax);
     
     % SaveData 
        filename       = FileName{n};
        SavePath{n}    = [filename num2str(labindex)];
        psave(filename,'F0','T0','f','t','p','dT','dF','-v6');                
        catch
         Key(n) = 0;
         Error{n} = err.message;
         disp(err.message)
        end       
    end
end

    data.detdata.EnergySignal   = EnergySignal;
    data.detdata.EnergySignaldB = EnergySignaldB;
    data.detdata.TimeSignal     = TimeSignal;
    data.detdata.TimeStep       = TimeStep;
    data.detdata.ClippingSignal = ClippingSignal;
    
    data.audioinfo.Key                  = Key;
    data.detdata.SavePath               = SavePath;
    data.detdata.DetectionTimeLength    = DetectionTimeLength;    
    data.detdata.Error                  = Error;
    
% 
%      E = sum(P(iF,:),1);        
%     %Clipping test
%         maxcount =10;
%         K = isclipped(audio,T,w,fs,maxcount); 
%     %
%         EnergySignal{n}    = E;
%         EnergySignaldB{n}  = 10*log10(E);
%         TimeSignal{n}      = T;
%         TimeStep(n)        = T(2)- T(1);
%         ClippingSignal{n}  = K;