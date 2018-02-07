function Data = localfeatext (Data,inFE,pkey)
    if nargin ==2
        pkey =0;
    end
    N = length(Data.AudioFileIndex);
    Data.FeatExt = inFE;
    % Feature Extraction data preparation
       info = featextparameters(inFE,Data);
    % Features Name    
    Data.FeatListStr =  Lfeatures()';
    Data.FeatList    = 1:size(Data.FeatListStr,2);
    
    % Feature Extraction parameters
        Data.TimeLim2 = info.det.T0;
        Data.SampLim2 = info.det.S0;
    if pkey
        ind = randperm(N);
    else
        ind = 1:N;
    end
    
    %Audio parameters
    path =  Data.AudioFilePath;
    fs   =  Data.Fs;
    % Audio segment extraction
    samples = info.det.S0;
    s0      = samples(:,1);
    frame   = info.stft;
    
    % Parameters
    D       = length(Lfeatures());
    eft     = cell(N,3);
    eftT0   = nan (N,1);
    X       = nan (N,D);
    F0      = nan (N,1);
    T0      = nan (N,1);
    theta   = -inFE.Threshold;
    snr     = nan (N,1);
    q       = nan (N,2);
    
    
    parfor n = ind
    disp(['Extracting features  ' num2str(n) ' of ' num2str(N)]);
    % Audio segment extraction
        audio   = wavread(path{n},samples(n,:));        
    % Time Frequency Representation
        [P,F,T,i0,j0,t0,~,~,~] = TFR(audio,frame,fs(n),n);    
    % Filtering
        S = TFRfilt(P,i0,j0,frame);
    %EFT Contour
        [e,f,t]     = EFTcontour (P,F,T,S);
        eft(n,:)    = {e,f,t};
        eftT0(n)    = t0 + s0(n)/fs(n);                
    % Feature Extraction
        X(n,:) = Lfeatures([],P,S,e,f,t);
        F0(n)  = F(i0);
        T0(n)  = T(j0);
    % Quality of the call
        snr(n) = SNR(P,F,S,theta);
        q(n,:) = Qmoments(P,F,[1 4])';
        
%         if pkey
%         
%                     E0 = sum(P(:));
%                     Q = 1-1/mean(kurtosis(10*log10(P),1));
% 
%                     figure(5)
%                     sigma = std((F*ones(1,size(P,2))).*P,1);
%                     sigma_max = (sum(P,1)*std(F));
%                     q1 = 1-(sigma_max-sigma)./sigma_max;
% 
%                     kappa = kurtosis(P,1);
%                     q2    =  1 - 1./kappa;
% 
%                     E = sum(P,1)/E0;
% 
%                     subplot (3,1,1)
%                     plot(q1);
%                     ylabel('std')
%                     xlabel(['mean = ' num2str(mean(q1))]);
%                     subplot (3,1,2)
%                     plot(q2);
%                     ylabel('kurtosis')
%                     xlabel(['mean = ' num2str(mean(q2))]);
%                     subplot (3,1,3)
%                     plot(E.*q1.*q2);
%                     ylabel('energy x std x kurtosis'); 
%                     xlabel(['mean = ' num2str(mean(E.*q1.*q2))]);
% 
%                     figure(1);clf
%                     hold on        
%                     surf(T,F,mag2db(P),'edgecolor','none'); axis tight; view(0,90);title(n)
%                     plot3(T(j0),F(i0),mag2db(P(i0,j0)),'.w','MarkerSize',20);
%                     x = info.stft.tplim(n,1); y = info.stft.fplim(n,1);
%                     dx = diff(info.stft.tplim(n,:)); dy = diff(info.stft.fplim(n,:));
%                     rectangle('position',[x, y, dx,dy],'EdgeColor','w','Linewidth',2);
%                     zlim([-1e3 0])
%                     xlabel(nnz(S)/nzmax(S)*info.stft.aframe)
%                     %title({'Energy frame = ' num2str(E0,'%2.5f'); 'Q = '  num2str(Q,'%2.5f')})
%                     title(['SNR = ' num2str(Data.SNR(n),'%2.5f') 'Q4 = '  num2str(Data.Q(n,2),'%2.5f') ' Q1 = '  num2str(Data.Q(n,1),'%2.5f')])
%                     
%                     hold off
% 
%                     figure(2);clf
%                     S0 = P*0; iP = P>db2pow(pow2db(P(i0,j0))- info.stft.thetha);
%                     S0(iP) = 1;
%                     image(T,F,S0*255); set(gca,'YDir','normal')
%                     Ef = sum(sum(S0.*P,1));
%                     Efstr = num2str(Ef/E0*100,'%%%1.2f');
%                     A     = nnz(S0)/nzmax(S0);
%                     Astr  = num2str(A*100,'%%%1.2f');
%                     title(['Energy - ' Efstr '    Area - ' Astr ]);        
% 
%                     figure(3);clf
%                     image(T,F,S*255); set(gca,'YDir','normal')
%                     Ef = sum(sum(S.*P,1));
%                     Efstr = num2str(Ef/E0*100,'%%%1.2f');
%                     A     = nnz(S)/nzmax(S);
%                     Astr  = num2str(A*100,'%%%1.2f');
%                     title(['Energy - ' Efstr ' Area - ' Astr ]);   
% 
%                     figure(4);clf
%                     hold on
%                     surf(T,F,mag2db(P.*S),'edgecolor','none'); axis tight; view(0,90);title(n)
%                     cmp = jet;
%                     set(gca,'color',cmp(1,:));
%                     rectangle('position',[x, y, dx,dy],'EdgeColor','w','Linewidth',2);
%                     plot3(t,f,mag2db(e),'-g','Linewidth',2);
%                     plot3(T(j0),F(i0),mag2db(P(i0,j0)),'.w','MarkerSize',20);
%                    zlim([-1e3 0]);
%                     hold off
%                     waitforbuttonpress
%         end
    end
    %EFT Contour
        Data.EFTcurve       = eft;
        Data.EFTt0          = eftT0;    
    % Feature Extraction
        Data.Data   = X;
        Data.F0     = F0;
        Data.T0     = T0;
    % Quality of the call
        Data.SNR    = snr;
        Data.Q      = q;    
end

    
