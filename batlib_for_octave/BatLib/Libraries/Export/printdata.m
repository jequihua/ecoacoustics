function printdata (data,objPrint,in,mode)


% Routine parameters
    N        = data.audioinfo.FilesAmount;
    Key      = data.audioinfo.Key;
    if isfield (data,'detdata')
        DetKey   = data.detdata.Key;
    end
    if isfield (data,'featext')
        FeatKey  = data.featext.Key;
    end
    
    if nargin ==3 || isempty(mode)
        mode = 'audiofile';
    end
        

%Page parameters
 
        %Plot A
        if isempty(in);
            in = 1:N;
        end
        
    % Audio file info     
    System   = data.path.Sytem;
    WinHeader= data.path.WinHeader;
    LinHeader= data.path.LinHeader;
    Path     = cell(N,1); 
    for i = in;
        switch mode
                case 'audiofile'
                    Path{i}     = changepath(data.audioinfo.Path{i},data.path);
                case 'calls'
                if isfield (data,'featext')
                    if FeatKey(i)
                    Path{i}    = changepath(data.featext.SavePath{i},data.path);
                    end
                end
        end
    end
    PrintFolder = changepath(data.path.Printouts,data.path);
    SaveFolder = fullfile(PrintFolder,[data.Title ' - ' mode]);
    SavePrefixName = fullfile(SaveFolder,[data.Title ' P']);
    if ~isdir(SaveFolder)
        mkdir(SaveFolder)
    end
    
    
    FileName = data.audioinfo.Name;
 
    % Audio signal parameters
    
    
    % Number of pages calculation
    
%     Ta = data.audioinfo.TimeLength(Key);
%     Ls = sum(ceil(Ta/tL(1)));
%     Tc = cellfun(@(x) sum(x),data.detdata.VADmethod.DetectionTimeLength);
%     Lc = sum(ceil(Tc(Key)/tL(2)));
%     
%     P  = ceil((Lc+Ls)/P.Nline);    
    
    
    l = 1;
    L = objPrint.trackNtotal;
    
    objPrint.audiofileNtotal = N;
    objPrint.pageNumber = 1;
    if isempty(objPrint.trackFsize)
                                autoFsize = true;
    end

    for n = in
        if Key(n)
            try
            sp = getspecparam (data,n,mode);
            skey  = data.detdata.Selection.Key{n}; 
            [P,F,T,x] = getspec(Path,sp,n,mode,skey);     
            M =  ceil(max(T)/objPrint.trackTsize);            
            
            switch mode
                        case 'audiofile'
                           if  autoFsize
                                if DetKey(n)
                                dt = getdetvar(data,n);
                                df  = 25e3;
                                maxF = max(ceil(max(dt.Ybox(:,2))/df)*df,50e3);
                                minF = 0;
                                objPrint.trackFsize = maxF;
                                objPrint.trackFlim  = [minF maxF];
                                else
                                    objPrint.trackFsize = data.audioinfo.SamplingFrequency(n)/2;
                                end
                            end
                        case 'calls'
                            if  autoFsize
                                df  = 25e3;
                                maxF = ceil (max(x.Fp)/df)*df;
                                minF = floor(min(x.Fp)/df)*df;
                                objPrint.trackFsize = maxF-minF;
                                objPrint.trackFlim  = [minF maxF];
                            end
                    end
            %objPrint.trackFsize = min(sp.fs/2,objPrint.trackFsize);
            objPrint.audiofileName = FileName{n};
            objPrint.audiofileNumber = n;
            objPrint.pageNote1 = sp;
            
            
            
            for m = 1:M               
                if l==1
                    h = objPrint.pageNumber;
                    objPrint.newfig(h);
                    objPrint.writenotes(h);
                end
                    
                    
                  objPrint.plottrack(P,F,T,l,m)                    
                    switch mode
                        case 'audiofile'
                            if DetKey(n)
                            dt = getdetvar(data,n);
                            c = objPrint.boxcolor;
                            objPrint.plotbox(dt.Xbox,dt.Ybox,c,l,m);
                            objPrint.plottext(dt.Xtx,dt.Ytx,dt.Fptx,c,l,m);
                            objPrint.plottext(dt.Xtx,dt.Ytx-objPrint.trackFsize/12,dt.dTtx,c,l,m);
                            objPrint.plotpoints(dt.Xp,dt.Yp,c,l,m);                    
                            if isfield(data,'classinfo')
                                    clname = 'RF_R1_Sp_F';
                                    cl = getclassdata(data,n,clname);
                                    objPrint.plottext(dt.Xtx,dt.Ytx*0+objPrint.trackFsize/12*2,dt.Ntx,c,l,m);  
                                    objPrint.plottext(dt.Xtx,dt.Ytx*0+objPrint.trackFsize/12,cl.Gatx,c,l,m);                                 
                            end

                            end
                        case 'calls'
                            if FeatKey(n)
                                
                                ft = getfeatvar(data,n);
                                ctx = [1 1 1]*0.95;
                                c = objPrint.boxcolor;
                                objPrint.plottext(x.T1+x.dT/2,minF*ones(size(x.T1))+ (maxF-minF)/8*2  ,ft.T0tx,c,l,m);
                                objPrint.plottext(x.T1+x.dT/2,minF*ones(size(x.T1))+ (maxF-minF)/8*3,ft.Fptx,c,l,m);
                                %objPrint.plottext(x.T1+x.dT/2,ft.Fbp+df/4,ft.Fptx,c,l,m);
                                %objPrint.plottext(x.T1+x.dT/2,x.F0-15e3,num2str(x.A,'%2.0f'),'w',l,m);
                                
                                %objPrint.plotpoints(x.Tp,x.Fp,c,l,m);
                                objPrint.plotpoints(ft.Tp+x.T0,ft.Fp,c,l,m);
                                for k = 1:length(ft.f)
                                    objPrint.plotpoints(ft.t{k}+x.T0(k),ft.f{k},'-',l,m);
                                end
                               
                                if isfield(data,'classinfo')
                                    clname = 'RF_R1_Sp_F';
                                    cl = getclassdata(data,n,clname);
                                    objPrint.plottext(x.T1+x.dT/2,minF*ones(size(x.T1))+ (maxF-minF)/8*1,cl.Gtx,c,l,m);                                    
                                end
                            end
                    end
                  
                if m==1 || l==1;
                   objPrint.writetitle(h,l); 
                end   
                if l==L || (m==M && n == in(end)) 
                    p = objPrint.pageNumber;
                    objPrint.savename = [SavePrefixName num2str(p,'%05.0f') ' - F' num2str(n,'%05.0f')];
                    objPrint.printpdf(h);
                    disp(objPrint.savename);
                    if objPrint.savekey
                    close(h);
                    end
                    objPrint.pageNumber = p+1; 
                    l = 1;
                else
                    l = l+1;
                end
            end
           
%             if DetKey(n)
%                 if FeatKey(n) 
%                     Lc = sum(ceil(Tc/tL(2)));
%                     for lc = 1:Lc
%                         
%                         
%                         
%                         if l==1
%                     
%                         elseif l==L
%                             l = 1;
%                             p = p+1;
%                         else
%                             l = l+1;
%                         end                
%                     end                    
%                 end                
%             end
            catch err
                disp(['Error in file:' num2str(n)]);
                disp(err.getReport)
            end
        end
    end
    
    
    
    function sp = getspecparam (data,n,mode)
        
        switch mode
            case 'audiofile' 
            sp.fs   = data.audioinfo.SamplingFrequency(n);
            if sp.fs<50e3
                sp.fs = sp.fs*10;
            end
            sp.wn   = data.detdata.WindowSize(n);
            sp.wt   = data.detdata.WindowTime(n);
            sp.wf   = data.detdata.WindowFunction;
            sp.ov   = round(data.detdata.Overlap*sp.wn); 
            sp.df   = sp.fs/sp.wn;
            sp.dt   = (sp.wn-sp.ov)/sp.fs;
            case 'calls'
            sp.fs   = data.audioinfo.SamplingFrequency(n);    
            sp.wn   = data.featext.STFT.WindowSize(n);
            sp.wt   = data.featext.STFT.WindowTime(n);
            sp.wf   = data.featext.STFT.WindowType;
            sp.ov   = round(data.featext.STFT.Overlap*sp.wn); 
            sp.df   = data.featext.STFT.FreqStep(n);
            sp.dt   = data.featext.STFT.TimeStep(n);
        end
            
        function dt = getdetvar(data,n)
                   k       = data.detdata.Selection.Key{n};
                   
                   if data.detdata.MethodKey == 2
                       method = 'VADmethod';
                   elseif data.detdata.MethodKey == 1
                       method = 'Emethod';
                   end
                   
                dt.Ybox    = data.detdata.(method).DetectionFreqInt{n}(k,:);
                dt.Xbox    = data.detdata.(method).DetectionTimeInt{n}(k,:);
                dt.Fptx     = num2str(data.detdata.(method).DetectionFreqPoint{n}(k,:)/1e3,'%2.1f');
                %dt.tx      = num2str(data.detdata.(method).FilledArea{n}(k,:),'%2.1f');
                dt.Ntx      = num2str(data.detdata.Selection.DetIndex{n}'); 
                dt.dTtx     = num2str(data.detdata.(method).DetectionTimeLength{n}(k,:)*1e3,'%2.1f');
                dt.Xtx     = data.detdata.(method).DetectionTimePoint{n}(k,:);
                dt.Ytx     = data.detdata.(method).DetectionFreqInt{n}(k,1);
                dt.Xp      = data.detdata.(method).DetectionTimePoint{n}(k,:);
                dt.Yp      = data.detdata.(method).DetectionFreqPoint{n}(k,:);
            function ft = getfeatvar(data,n)
                k       = data.detdata.Selection.Key{n};
                ft.f = data.featext.EFT.Curve{n}(k,2);
                ft.t = data.featext.EFT.Curve{n}(k,3);
                
                X = data.featext.EFT.Features{n}(k,:);
                        
                ft.Fk = X(:,strcmp(data.featext.EFT.FeatNames,'Frequency of the knee')');
                ft.Tk = X(:,strcmp(data.featext.EFT.FeatNames,'Time of the knee')');
                
                ft.Fc = X(:,strcmp(data.featext.EFT.FeatNames,'Characteristic Frequency')');
                ft.Tc = X(:,strcmp(data.featext.EFT.FeatNames,'Characteristic Frequency Time')');
                
                ft.Fb = X(:,strcmp(data.featext.EFT.FeatNames,'Bottom Frequency')');
                ft.Ftop = X(:,strcmp(data.featext.EFT.FeatNames,'Top Frequency')');
                                
                ft.Fp   = X(:,strcmp(data.featext.EFT.FeatNames,'Peak Frequency')');
                ft.Tp   = zeros(size(X,1),1);
                for i = 1:size(X,1)
                ft.Tp(i)   = ft.t{i}(find(ft.f{i} == ft.Fp(i),1,'first'));
                end
                
                ft.T0tx =num2str(data.featext.EFT.T0{n}(k,:),'%2.2f');
                ft.Fptx =num2str(X(:,strcmp(data.featext.EFT.FeatNames,'Peak Frequency')')/1e3,'%2.1f');
                
                function cl = getclassdata(data,n,clname)
                     %k       = data.detdata.Selection.Key{n};
                     K       = data.detdata.Selection.Amount(n);
                     %ik      = data.detdata.Selection.DetIndex{n};
                     Gtx      = data.classinfo.(clname).Prediction{n};
                     cl.Gtx   = Gtx;
                     cl.Gatx  = cell(K,1);
                     for i = 1:K
                         itx =  regexp(Gtx{i},'([aeiou][^aeiou]+)','end');
                         cl.Gatx{i}  = Gtx{i}(1:itx(2));
                     end
                     
            function [P,F,T,x] = getspec(Path,sp,n,mode,skey)
                
                switch mode
                    case 'audiofile'                      
                      audio   = audioread (Path{n});
                      [P,F,T] = spec(audio,sp.wn,sp.ov,[],sp.fs,sp.wf,'tftb');
                      x      =0;
                    case 'calls'
                        load(Path{n});
                        C = selectC(C,skey);
                        [P,F,T,x] = callspec(C);
                end
                
                function C = selectC(C,k)
                    K0    = C.NumObjects;
                    K     = nnz(k);
                    cname = fieldnames(C);
                    for i = 1:size(cname,1);
                        if size(C.(cname{i}),1) == K0
                            C.(cname{i}) = C.(cname{i})(k,:);
                        end
                    end
                    C.NumObjects = K;
                    
            function [P,F,T,r] = callspec(C)
                
                
                x1   = C.BoundingBox(:,1)+ 1/2;
                x2   = C.BoundingBox(:,3)+ x1-1;
                dx   = C.BoundingBox(:,3);
                
                 a   = round(10e-3/C.dT);
                dxa  = ceil(dx/a)*a;
                %dxa(dx<a/2)= round(a/2);                
                j0   = cumsum([1;dxa(1:end-1)])+round((dxa-dx)/2);
                
                i0   = C.F0/C.dF+1;
                %y1   = C.BoundingBox(:,2)+ 1/2+I0;
                %y2   = C.BoundingBox(:,4)+ y1-1;
                %dy   = y2-y1+1; 

                Imax = max(C.ImageSize(:,1)+i0);
                Jmax = sum(dxa);                
                minP = 10^floor(log10(min(vertcat(C.PixelValues{:}))/4));
                P    = ones(Imax,Jmax)*minP;
            for k =1:C.NumObjects                
                i   = C.ImageSize(k,1);
                j   = C.ImageSize(k,2);
                p0  = ones(i,j)*minP;   
                p0(C.PixelIdxList{k})= C.PixelValues{k};
                p   = p0(:,x1(k):x2(k));
                P (i0(k):(i0(k)+i-1), j0(k):(j0(k)+dx(k)-1))= p;
            end
            
                F  = (1:Imax)'*C.dF;
                T  = (0:Jmax-1)*C.dT;
               
              ri   = zeros(C.NumObjects,1);
              rj   = zeros(C.NumObjects,1);
              for k =1:C.NumObjects   
                [Ri,Rj] = ind2sub([i j],C.RefPoint(k));
                ri(k)   = (Ri+i0(k)-1);
                rj(k)   = (Rj+j0(k)-x1(k)-1);            
              end
              
              r.Tp = rj*C.dT;
              r.Fp = ri*C.dF;
              r.T0 = (j0-x1)*C.dT-C.T0;
              r.F0 = i0*C.dF;
              r.T1 = (j0)*C.dT;              
              r.dT = (dx)*C.dT;
              r.dF = (C.BoundingBox(:,2)+ C.BoundingBox(:,4))*C.dF;              
              r.A  = C.Area;
              
            
            
    
            