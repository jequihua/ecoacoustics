function data = featext(data)
        
    % Routine parameters
    N        = data.audioinfo.FilesAmount;
    Key      = data.detdata.Key;
    Error    = data.detdata.Error;
    %SelAmn   = data.detdata.Selection.Amount;
    switch data.detdata.MethodKey
        case 1
            method = 'Emethod'; 
        case 2
            method = 'VADmethod';
    end
    Ndet   = data.detdata.(method).NumberOfDetections;
    
    
    % Audio file info    
    
    Name     = data.audioinfo.Name;
    Path     = changepath(data.audioinfo.Path,data.path);
    SavePath = cell(N,1);
    SaveMode = data.featext.SaveMode;
    for i = 1:N
        if Key(i)
            FileTag  = ['F' num2str(data.audioinfo.FileIndex(i),'%06d') ' - ' Name{i}];            
            filepath = fullfile(changepath(data.path.Sonograms,data.path),FileTag);
            SavePath{i} = filepath;
        else
            SavePath{i} = '';
        end
    end 
    
    try
        Pool     = gcp('nocreate');
        if isempty(Pool)
        poolsize = 1;
        else
        poolsize = Pool.NumWorkers;
        end
    catch
        poolsize=12;
    end
    
    DET  = cell(N,1);
    %PSDP = cell(N,1);
    STFT  = cell(N,1);       
    
    
    WindowTime = zeros(N,1);
    WindowSize = zeros(N,1);
    TimeStep   = zeros(N,1);
    FreqStep   = zeros(N,1);
    
    for n = 1:N;
        if Key(n)
        %Detection parameters    
        det  = detparameters(data,n);
        %PSD parameters
        %psdp  = psdparameters(data,psdp,Fs(n),det);
        %STFT parameters
        stft = stftparameters(data,n);
        %Selection of detection
        m  = true(Ndet(n),1);
        %m    = data.detdata.Selection.Key{n};
        % Samples limits
        [DET{n},STFT{n}] = getsamplelimits(det,stft,m);
        
        WindowTime(n) = stft.wn/stft.fs;
        WindowSize(n) = stft.wn;
        FreqStep(n)   = stft.fstep;
        TimeStep(n)   = stft.tstep;
        end
    end

    EFTcurve   = cell(N,1);
    EFTt0      = cell(N,1);
    EFTfeat    = cell(N,1);
    EFTfeatnames   = Lfeatures;
    D = size(EFTfeatnames,1);

    mtotal = sum(Ndet);
    itime = now;
    hxs = 1/(60*60*24);
    
parfor n = 1:N;
    if Key(n)
        path    = Path{n};    
        det     = DET{n};
        stft    = STFT{n};

        samples = det.S0;
        T0      = det.T0;    
        M       = Ndet(n); 
        eft     = cell(M,3);
        eftT0   = zeros(M,1);
        X       = nan(M,D);
        if ~isempty(SaveMode)
            CP.dF       = stft.fstep;
            CP.dT       = stft.tstep;
            CP.Connectivity = stft.conn;
            CP.NumObjects = M;    
        end

        disp(['Feature extraction. File ' num2str(n) ': ' Name{n} ' - ' num2str(M) ' Calls' ])
        tic
        try
            for m = 1:M;
                try
                %disp(['EFT Feature extraction. File ' num2str(n) ': Call ' num2str(m) ' of ' num2str(M)])

                %Audio reading            
                audio   = audioread(path,samples(m,:));               
    %               %STFT features
                %Spectrogram frame
                [P,F,T,i0,j0,~,~,~] = stftframe(audio,stft,m);

                %Get signal distribution
                [S,cp] = getsignaldist(P,i0,j0,stft,1);

                if ~isempty(SaveMode)
                    cp.T0 = T(1);
                    cp.F0 = F(1);                
                    CP = packcp(CP,cp,m,M);
                end            
                %EFT Contour
                [e,f,t,c]     = EFTcontour2 (P,F,T,S,0.5e-3);
                if length(c.Tc)< 4
                    disp('EFTcontour2 adjusted')                    
                end
                
                eft(m,:)    = {e,f,t};
                eftT0(m,1)  = T(1) + T0(m);
                X(m,:) = Lfeatures ([],P,S,e,f,t,c);
                %plotcall(P,F,T,S,Pi,Fi,Ti,i0,j0,e,f,t,stft,m)
    %                 if any(isnan(e))|| any(isnan(f))|| any(isnan(t))
    %                     plotcall(P,F,T,S,Pi,Fi,Ti,i0,j0,e,f,t,stft,m)
    %                     disp(find(isnan(e)));
    %                 end
                catch err
                    disp([n m M]);
                    disp(err.getReport)                    
                end
            end

            % Disp Info
                    if ~isempty(SaveMode)                   
                        savefeatvar(P,F,T,CP,SavePath{n},SaveMode)
                    end

                    mtime = toc/M;               
                    ftime =  datestr(datenum(itime)+datenum(mtime*mtotal*toc*hxs));
                    dispstr = ['End time estimation '  ftime];
                    disp(dispstr)

                    EFTcurve{n} = eft;
                    EFTt0{n}    = eftT0;
                    EFTfeat{n}  = X;

        catch err
                    Key(n) = 0;
                    Error{n} = err;
                    disp(err.getReport)
        end
    end
end
       
        
        data.featext.Error                  = Error;
        data.featext.Key                    = Key;
        try
        data.featext.SavePath               = changepath(SavePath,data.path,computer);        
        end
        data.featext.STFT.WindowTime =   WindowTime;
        data.featext.STFT.WindowSize =   WindowSize;
        data.featext.STFT.TimeStep   =   TimeStep;
        data.featext.STFT.FreqStep   =   FreqStep;
        
        data.featext.EFT.Curve              = EFTcurve;
        data.featext.EFT.T0                 = EFTt0;
        data.featext.EFT.Features           = EFTfeat;
        data.featext.EFT.Nfeat              = D;
        data.featext.EFT.FeatNames          = EFTfeatnames;
        
        function savefeatvar(p,f,t,C,filename,mode)
switch mode
    case 's'
        save(filename,'f','t','p','-v6');
    case 'd'
        save(filename,'C','-v6'); 
    case 's+d'
        save(filename,'f','t','p','C','-v6');
    otherwise            
end
function CP = packcp(CP,cp,m,M)    
    
    if m == 1        
        CP.ImageSize        = nan(M,2);
        CP.PixelIdxList0    = cell(M,1); 
        CP.PixelIdxList     = cell(M,1);
        CP.Area0            = nan(M,1);
        CP.Area             = nan(M,1);
        CP.BoundingBox      = nan(M,4);
        CP.PixelValues      = cell(M,1);
        CP.RefPoint         = nan(M,1);
        CP.T0               = nan(M,1);
        CP.F0               = nan(M,1);
    end
    
        CP.ImageSize(m,:)        = cp.ImageSize;
        CP.PixelIdxList0(m)      = cp.PixelIdxList0;
        CP.PixelIdxList(m)       = cp.PixelIdxList;
        CP.Area0(m,:)            = cp.Area0*CP.dF*CP.dT;
        CP.Area(m,:)             = cp.Area*CP.dF*CP.dT;
        CP.BoundingBox(m,:)      = cp.BoundingBox;
        CP.PixelValues{m}        = cp.PixelValues;
        CP.RefPoint(m,:)         = cp.RefPoint;
        CP.T0(m,:)               = cp.T0;
        CP.F0 (m,:)              = cp.F0;
    
    
    

       