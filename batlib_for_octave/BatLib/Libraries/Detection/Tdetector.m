function [D,onset,offset,Ymax,iYmax] = Tdetector(Y,theta)
            
        D       = Y >= theta;
        onset   = [0  diff(D)==1];
        offset  = [diff(D)==-1 0];

        if sum(onset) < sum(offset)
            onset (1) = 1;
        elseif sum(onset) > sum(offset)
            offset (end) = 1;
        end
        Onset  = find(onset)';
        Offset = find(offset)';
        
        P = length(Onset);
        iYmax = zeros(P,1);
        Ymax  = nan(P,1);
        for p = 1:P;
            [ymax iymax] = max(Y(Onset(p):Offset(p)));
            iYmax(p) = iymax + Onset(p)-1;
             Ymax(p) = ymax; 
        end