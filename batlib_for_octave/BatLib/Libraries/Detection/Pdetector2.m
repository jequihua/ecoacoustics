function [D, Onset, Offset,Ymax,iYmax] = Pdetector2(Y,alpha,delta)
    
            [peaks,valleys] = peakdet(Y,alpha);            
            
            iYmax = peaks(:,1);     
            Ymax  = peaks(:,2);
            
            N =  length(Y);
            P =  length(Ymax);
            
            [Ymin0 iYmin0] =  min(Y(1:iYmax(1)));
            iYmin = [iYmin0; valleys(:,1)]; 
            %Ymin  = [ Ymin0; valleys(:,2)];
            
            if size(iYmax,1) == size(iYmin,1)
                   [yminf iyminf] = min (Y(iYmax(P):N));
                   iYmin = [iYmin; iyminf + iYmax(P)- 1];
                   %Ymin  = [Ymin; yminf];
            end
            
            
              
            Onset  = max(iYmax - floor(delta/2),1);
            Offset = min(iYmax + floor(delta/2),N);
            
               
            deltaYmax = abs(diff([1; iYmax; N]));
            i = find(deltaYmax < delta);
            I = length(i);
            
            Onset (i(1:I-1)) = iYmin(i(1:I-1));
            Offset(i(2:I)-1) = iYmin(i(2:I));
            
            D = zeros(1,N);
            for p = 1:P;                
                D (Onset(p):Offset(p))= 1;                             
            end