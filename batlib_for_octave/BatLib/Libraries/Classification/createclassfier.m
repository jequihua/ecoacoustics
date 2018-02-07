function data = createclassfier (data)

    
    Key      = data.audioinfo.Key;
    Name     = data.audioinfo.Name;
    SavePath = data.detdata.SavePath;
    System   = data.path.Sytem;
    WinHeader= data.path.WinHeader;
    LinHeader= data.path.LinHeader;    
    Error    = data.detdata.Error;
    
    iN = ~(data.detdata.Selection.Amount==0);
    
    for c = 1: length(data.classinfo.RF)
        disp(data.classinfo.RF(c).Name);
        R = data.classinfo.RF(c).Region;
        ir = find(sum(data.classinfo.Region(:,find(R)),2)>0);
        iR = false(size(iN));
        for n = 1:data.audioinfo.FilesAmount
           if any(data.classinfo.SpeciesIndex(n)==ir);
               iR(n)=1;
           end
        end
        
        M = data.classinfo.RF(c).RecordingMode;
        if ~isempty(M)
            iM = false(size(iN));
            for n = 1:data.audioinfo.FilesAmount
               if strcmp(data.classinfo.RecordingMode{n},M);
                   iM(n)=1;
               end
            end
        else
            iM = iN;
        end

        iK = iN&iR&iM;
    
        N = sum(data.detdata.Selection.Amount(iK));
        D = data.featext.EFT.Nfeat;
        Y = cell(N,1);
        X = nan(N,D);
        n0 = 1;
        for i = find(iK)'
            M = data.detdata.Selection.Amount(i);
            Y(n0:n0+M-1) = data.classinfo.SpeciesStr(i);
            X(n0:n0+M-1,1:D) = data.featext.EFT.Features{i};
            n0 = n0+M;
        end
        opts = statset('UseParallel',true,'Display','iter');
        T = TreeBagger(500,X,Y,'Method','classification','Options',opts,'Prior','Uniform','NVarToSample',7);
        data.classinfo.RF(c).class = T.compact;        
    end
    