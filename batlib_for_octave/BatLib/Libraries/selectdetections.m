function data = selectdetections(data,nfile)

     minN = data.detdata.Selection.AmountLim(1);
     maxN = data.detdata.Selection.AmountLim(2);
     
     SelMode = data.detdata.Selection.Mode; 
     
     switch data.detdata.MethodKey
         case 1
             detdata = data.detdata.Emethod;
         case 2
             detdata = data.detdata.VADmethod;             
     end
     
        mindTD  = data.detdata.Selection.CallDurLim(1);
        maxdTD  = data.detdata.Selection.CallDurLim(2);
        
        minT    = data.detdata.Selection.TimeLim(1);
        maxT    = data.audioinfo.TimeLength - data.detdata.Selection.TimeLim(2);
        if isfield(data.detdata.Selection,'DetectionLRLim')
            minLR =  data.detdata.Selection.DetectionLRLim(1);
            maxLR =  data.detdata.Selection.DetectionLRLim(2);
        end
     switch SelMode
            case 'MaxE'
                
            case 'Area'                
              minA    = data.detdata.Selection.AreaLim(1);
              maxA    = data.detdata.Selection.AreaLim(2);
            case 'FreqFilter'
                minA    = data.detdata.Selection.AreaLim(1);
                maxA    = data.detdata.Selection.AreaLim(2);
                
            case 'BestSeq'
              minM    = data.detdata.Selection.GroupSize(1);
              maxM    = data.detdata.Selection.GroupSize(2);
                                          
              mindTP  = data.detdata.Selection.CallInterval(1);
              maxdTP  = data.detdata.Selection.CallInterval(2);
                            
              minG    = data.detdata.Selection.GammaFactor(1);
              maxG    = data.detdata.Selection.GammaFactor(2);
              Gmean   = data.detdata.Selection.GamaFactorMean;
             
     end
     
     
     
        DKey = data.detdata.Key;
         try 
            FKey = data.featext.Key;
         catch
            FKey = DKey;
         end
         
    N        = data.audioinfo.FilesAmount;
    if nargin <2
        k = find(DKey & FKey)';        
        SelKey   = cell(N,1);
        SubIndex = cell(N,1);
        Index    = cell(N,1);
        Amount   = zeros(N,1);
    else
        k = nfile;
        SelKey   = data.detdata.Selection.Key;
        SubIndex = data.detdata.Selection.DetIndex;
        Index    = data.detdata.Selection.SelIndex;
        Amount   = data.detdata.Selection.Amount;
    end
    
    
            
     for n = k
        Ndet  = detdata.NumberOfDetections(n);
        if Ndet >0 
            if data.detdata.Selection.ClippingFilter;
                iC     = ~detdata.DetectionsClipSignal{n};
            else
                iC     = true(Ndet,1);
            end

            dTD     = detdata.DetectionTimeLength{n};
            idTD    = mindTD<=dTD & dTD <=maxdTD;

            T       = detdata.DetectionTimePoint{n};        
            iT      = minT<=T & T <=maxT(n);
            
            if data.detdata.Selection.FiltNaNFeatures
                try                                
                iX      = any(~isnan(data.featext.EFT.Features{n}),2);
                end
                if isempty(iX)
                    iX      = true(Ndet,1);
                end
            else
                iX      = true(Ndet,1);
            end
            
            if isfield(data.detdata.Selection,'DetectionLRLim')
                LR    = detdata.DetectionLRValue{n}(:,2);
                iLR   = minLR<=LR & LR <=maxLR;
            else
                iLR   = true(Ndet,1);
            end
            
            if isfield(data.detdata.Selection,'CricketFilter')
                if data.detdata.Selection.CricketFilter;
                    if isfield(data.detdata,'CricketFilter')
                        if data.detdata.CricketFilter.Key(n);
                             dT    = detdata.DetectionTimeInt{n};
                             dF    = detdata.DetectionFreqInt{n}; 
                             dTc   = data.detdata.CricketFilter.DetectionTimeInt{n};
                             dFc   = data.detdata.CricketFilter.DetectionFreqInt{n};
                             nC    = size(dFc,1);
                             iCF   = true(Ndet,1);
                            for i = 1:Ndet
                                for j = 1: nC
                                    if dTc(j,1) < dT(i,2) &  dT(i,1) < dTc(j,2) 
                                        if dFc(1) < dF(i,2) &  dF(i,1) < dFc(2)
                                            iCF(i) = 0;
                                        end
                                    end
                                end
                            end
                        else
                            iCF  = true(Ndet,1);
                        end
                    else
                           iCF   = true(Ndet,1);
                    end
                    disp(sum(iCF));
                else
                    iCF = true(Ndet,1);
                end
            end
 
            switch SelMode
                case 'MaxE'
                E   = detdata.DetectionPowPoint{n};
                selkey = false(Ndet,1);            
                    if Ndet >= minN
                        [~,iEs] = sort(E,'descend');
                        iTDC = find(idTD & iC & iT);         
                        D = iEs(find(ismember(iEs,iTDC),maxN,'first'));
                        selkey(D)  = 1;
                        %Serie   = Key;
                        %S       = 1;
                    end
                case 'Area'

                    A   = detdata.FilledArea{n}; 
                    iA  = minA<=A & A <=maxA;
                    
                    Fmin  = data.detdata.Selection.FrequencyFilter(1);
                    Fmax  = data.detdata.Selection.FrequencyFilter(2);
                    
                    %    iFp = strcmp(data.featext.EFT.FeatNames,'Peak Frequency');
                    %    Fp  = data.featext.EFT.Features{n}(:,iFp);
                    Fp    = data.detdata.VADmethod.DetectionFreqPoint{n};
                    iF    = Fmin <= Fp & Fp <= Fmax; 
                    
                    selkey = idTD & iC & iT & iA & iX & iLR & iF & iCF;
                    %Serie   = Key;
                    %S       = 1;
                case 'FreqFilter'
                    
                    A   = detdata.FilledArea{n}; 
                    iA  = minA<=A & A <=maxA;
                    %selkey = idTD & iC & iT & iA;
                    
                    Fmin = data.detdata.Selection.FrequencyFilter(n,1);
                    Fmax = data.detdata.Selection.FrequencyFilter(n,2);
                    Fp    = data.detdata.VADmethod.DetectionFreqPoint{n};
                    iF    = Fmin <= Fp & Fp <= Fmax; 
                    selkey = idTD & iC & iT & iF & iA & iX;
                    
                case 'FreqFilter + maxE'
                    Fmin = data.detdata.Selection.FrequencyFilter(1);
                    Fmax = data.detdata.Selection.FrequencyFilter(2);
                    Fp    = data.detdata.Emethod.DetectionFreqPoint{n};
                    iF    = Fmin <= Fp & Fp <= Fmax;
                                        
                    iD    = idTD & iC & iT & iF & iX;
                    selkey = false(Ndet,1);
                    if sum(iD)> minN
                        E   = detdata.DetectionPowPoint{n};
                        [~,iEs] = sort(E,'descend');
                        D = iEs(find(ismember(iEs,find(iD)),maxN,'first'));
                        selkey(D) = 1;
                    else
                        selkey = iD;
                    end 
                case 'BestSeq'
                    Fp    = data.detdata.Emethod.DetectionFreqPoint{n};
                    mFp   = mean(Fp);
                    sFp   = std(Fp)/2;
                    Fmin  = max(data.detdata.Selection.FrequencyFilter(1),mFp - sFp);
                    Fmax  = min(data.detdata.Selection.FrequencyFilter(2),mFp + sFp);                    
                    iF    = Fmin <= Fp & Fp <= Fmax;
                    
                    iE      = data.detdata.Emethod.DetectionPoint{n}(:,2);
                    E       = data.detdata.Emethod.DetectionSignal{n}(iE);
                    D       = iC & idTD& iT & iF;
                    [M,M12,T12] = groupdet(T(),Gmean,minG,maxG,maxdTP);
                    iM = minM <= M & M <= maxM;
              
                      for i= 1:length(M)
                          j = M12(i,1):M12(i,2); 
                          if any(~D(j))
                              iM(i) = 0;
                          end                  
                      end             

                      m   = M(iM);
                      m12 = M12(iM,:);

                      selkey    = false(N,1);
                      Serie  = zeros(N,1);
                      %S      = 0;

                      if sum(m)> maxN
                            Emean = zeros(length(m),1);
                            for i = 1:length(m)
                                j = m12(i,1):m12(i,2);
                                Emean(i) =  mean(db2mag(E(j)));
                            end                            
                            [~,iEs] = sort(Emean,'descend');
                            ms   = m(iEs);
                            m12s = m12(iEs,:);

                            Ncum = cumsum(ms);
                            S = find(Ncum>=maxN,1,'first');
                            for i = 1: S-1
                                k = m12s(i,1):m12s(i,2);
                                selkey(k)  = 1;
                                Serie(k)= i; 
                            end
                                k = m12s(S,1):m12s(S,2)-(Ncum(S)-maxN);
                                selkey(k) = 1;
                                Serie(k)   = S;
                      elseif sum(m)>= minN
                            S = length(m);
                            for i = 1: S
                                k = m12(i,1):m12(i,2);
                                selkey(k)  = 1;
                                Serie(k)= i;
                            end                            
                      else
                          if Ndet >= minN
                            [~,iEs] = sort(E,'descend');
                            iTDC = find(idTD & iC & iT);         
                            D = iEs(find(ismember(iEs,iTDC),maxN,'first'));
                            selkey(D)  = 1;                        
                          else
                            selkey(D)  = 1;
                          end
                            Serie   = selkey;
                            S       = 1;                          
                          disp('No call group was found')
                      end
            end
           
            K = nnz(selkey);
            index = zeros(Ndet,1);
            index(selkey) = 1:K; 

            disp(['Selecting Calls. File ' num2str(n) ': ' data.audioinfo.Name{n} '(' num2str(K-Ndet) ' Calls)' ]) 
            SelKey{n}     = selkey;
            SubIndex{n}   = find(selkey');
            Amount(n)     = K;
            Index{n}      = index;
        end
     end
            data.detdata.Key                  = Amount > 0;
            data.detdata.Selection.Key        = SelKey;
            data.detdata.Selection.DetIndex   = SubIndex;
            data.detdata.Selection.Amount     = Amount;
            data.detdata.Selection.SelIndex   = Index;