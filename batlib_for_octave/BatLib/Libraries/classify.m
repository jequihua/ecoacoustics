function data = classify(data,obj)

    Error                   = data.detdata.Error;
    Path = changepath(data.audioinfo.Path,data.path);    
    Fs   = data.audioinfo.SamplingFrequency;
    N    = data.audioinfo.FilesAmount;
    Name     = data.audioinfo.Name;
    
    Pool                    = gcp('nocreate');
    if isempty(Pool)
    poolsize = 1;
    else
    poolsize = Pool.NumWorkers;
    end
    
    
    X = cell(nnz(data.detdata.Selection.Amount));
    switch obj.FeatType
        case 'detdata'
            for n = find(data.detdata.Selection.Amount)'
                X{n} = detfeat(data,n,obj);
            end
        case 'featdata'
            [~,jD,~]    = intersect(data.featext.EFT.FeatNames,obj.FeatNames);
            j           = false(1,data.featext.EFT.Nfeat);
            j(jD)       = 1;
            for n = find(data.detdata.Selection.Amount)'
                i    = data.detdata.Selection.Key{n};
                X{n} = data.featext.EFT.Features{n}(i,j);
            end
        case 'detdata+featdata'
            [~,jD,~]    = intersect(data.featext.EFT.FeatNames,obj.FeatNames);
            j           = false(1,data.featext.EFT.Nfeat);
            j(jD)       = 1;
            for n = find(data.detdata.Selection.Amount)'                
                i    = data.detdata.Selection.Key{n};
                X1   = data.featext.EFT.Features{n}(i,j);
                X2   = detfeat(data,n,obj);
                X{n} = [X1 X2];
            end            
        end
    
    
    Prediction = cell(N,1);
    PostProb  = cell(N,1);
    C = obj.class;
    key =  data.detdata.Selection.Amount;

    
    itime = now;
    hxs = 1/(60*60*24);
    mtotal = sum(data.detdata.Selection.Amount);
parfor n = 1:N
    tic
    if key(n)
        try
        disp(['Classification. File ' num2str(n) ': ' Name{n}])
        [prediction,postprob] = predict(C,X{n});
        postprob = max(postprob,[],2);
        Prediction{n} = prediction;
        PostProb {n}  = postprob;
            % Disp Info
                mtime = toc/key(n);               
                ftime =  datestr(datenum(itime)+datenum(mtime*mtotal*toc*hxs));
                dispstr = ['End time estimation '  ftime];
                disp(dispstr)
        catch err
        Error{n} = err.message;
        disp(err.message)
        end
    end
end
data.detdata.Error                      = Error;
data.classinfo.(obj.Name).Prediction    = Prediction;
data.classinfo.(obj.Name).PostProb      = PostProb;
data.classinfo.(obj.Name).Region        = obj.Region;
data.classinfo.(obj.Name).RecordingMode = obj.RecordingMode;

function x = detfeat(data,n,obj)
        D = data.detdata.VADmethod;
        J = length(obj.FeatNames);
        f = obj.FeatNames;
        i = data.detdata.Selection.Key{n};
        I = data.detdata.Selection.Amount(n);
        
        x = nan(I,J);
     for j = 1:J
        switch f{j}
            case 'DetectionTimeLength'
                x(:,j) = D.(f{j}){n}(i,:)*1e3;
                
            case 'DetectionFreqLength'
                x(:,j) = D.(f{j}){n}(i,:)*1e-3;
                
            case 'DetectionFreqPoint'
                x(:,j) = D.(f{j}){n}(i,:)*1e-3;
                
            case 'DetectionPowPoint'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'DetectionPowRng'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'DetectionPowMean'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'Eccentricity'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'Orientation'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'Solidity'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'MajorAxisLength'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'MinorAxisLength'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'FilledArea'
                x(:,j) = D.(f{j}){n}(i,:);
               
            case 'TimeCentroid';
                x(:,j) = D.Centroid{n}(i,1)*1e3;
                
            case 'FreqCentroid';
                x(:,j) = D.Centroid{n}(i,2)*1e-3;

            case 'TimeWeightedCentroid'
                x(:,j) = D.WeightedCentroid{n}(i,1)*1e3;
                
            case 'FreqWeightedCentroid'
                x(:,j) = D.WeightedCentroid{n}(i,2)*1e-3;

            case 'DetectionLRValue'
                x(:,j) = D.(f{j}){n}(i,2)
               
        end
     end