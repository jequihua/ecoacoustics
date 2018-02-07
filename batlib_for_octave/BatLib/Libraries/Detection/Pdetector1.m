function [D, Onset, Offset,Ymax,iYmax] = Pdetector1(Y,alpha,beta)
    
            [peaks,valleys] = peakdet(Y,alpha);
            
            if isempty(peaks)
            end
            
            iYmax = peaks(:,1);     
            Ymax  = peaks(:,2);
            
            N = length(Y);
            P = length(Ymax);
            
            if isempty(valleys)
                if size(Ymax,1)== 1;
                    iv = find(Y(iYmax:N)<= (Ymax-alpha),1,'first') +iYmax -1;
                    valleys(1,1) = iv;
                    valleys(1,2) = Y(iv);
                end
            end
            
            [Ymin0 iYmin0] =  min(Y(1:iYmax(1)));
            iYmin = [iYmin0; valleys(:,1)]; 
            Ymin  = [ Ymin0; valleys(:,2)];
            
            if size(iYmax,1) == size(iYmin,1)
                   [yminf iyminf] = min (Y(iYmax(P):N));
                   iYmin = [iYmin; iyminf + iYmax(P)- 1];
                   Ymin  = [Ymin; yminf];
            end
            
            
            % Determine High of every peak P
            y0 = max([Ymin(1:P) Ymin(2:P+1)],[],2);            
            dy = abs(Ymax-y0);
            
            theta = Ymax - beta*dy;
            
            Onset   = zeros(P,1); 
            Offset  = zeros(P,1);
            D       = zeros(1,N);
            
          for p = 1:P                
                Onset(p)  = find( Y(iYmin(p):iYmax(p))  <= theta(p),1,'last') +iYmin(p)-1;
                Offset(p) = find( Y(iYmax(p):iYmin(p+1))<= theta(p),1,'first')+iYmax(p)-1;
                D (Onset(p):Offset(p))= 1;                             
          end