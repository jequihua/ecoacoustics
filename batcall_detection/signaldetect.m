function out = signaldetect(Sx,Sxx,theta,t,f,prm)
        D       = Sx >= theta;
        
        onset   = [0  diff(D)==1];
        offset  = [diff(D)==-1 0];

        if sum(onset) < sum(offset)
            onset (1) = 1;
        elseif sum(onset) > sum(offset)
            offset (end) = 1;
        end
        
        Onset  = find(onset)';
        Offset = find(offset)';
        
        if prm.vad.ne == 0;
            if ~isempty(Onset)
                if Onset(1) < round(prm.ne.td/prm.vad.dt)
                D(Onset(1):Offset(1))=0;
                Onset(1) = [];
                Offset(1)= [];
                end
            end
        end
        
        P = length(Onset);
        if P>0
        iPmax = zeros(P,1);
        Pmax  = zeros(P,1);
        
            for p = 1:P;
                [pmax,ipmax] = max(Sx(Onset(p):Offset(p)));
                iPmax(p) = ipmax + Onset(p)-1;
                 Pmax(p) = pmax; 
            end
        
        iPx    = [Onset iPmax Offset];
        Px     = Pmax;
        
        [Py,iPy] = max(Sxx(:,iPmax),[],1);
        end

        iPt = t(iPx);
        out.NumberOfDetections    = P;
        out.DetectionTimeInt      = [iPt(:,1) iPt(:,3)];
        out.DetectionTimeLength   = iPt(:,3) - iPt(:,1);
        out.DetectionTimePoint    = iPt(:,2);
        
        iPf = f(iPy);
        out.DetectionFreqInt      = Py;
        out.DetectionFreqPoint    = iPf;
        
        out.DetectionSmpPoint     = [iPx(:,2) iPy'];
        out.DetectionPowPoint     = [Px Py'];
        out.DetectionSignal       = D;

        