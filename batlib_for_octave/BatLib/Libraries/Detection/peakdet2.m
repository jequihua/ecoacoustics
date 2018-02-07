function [Yp,iYp,hYp,wYp,dYp,y] = peakdet2(Y,hYpmin,nYp,sortYp,dYpType,wYpLevel)
%peak detectiion base on peakdet algotirthm with additional heigth 


if nargin <6
        wYpLevel = 0.5;
end
if nargin <5
        dYpType = 'height';
end
if nargin <4
        sortYp  = 'none';
end
if nargin <3
        nYp     = [];
end
if nargin <2
        hYpmin   = nanstd(Y)/2;
end

if isempty(hYpmin)
        hYpmin   = nanstd(Y)/2;
end




        N = length(Y);
        [pks,vlls] = peakdet(Y,hYpmin);
        
        iYp = pks(:,1);     
         Yp = pks(:,2);
         P  = length(Yp); 
         
         if isempty(vlls)
                if P == 1;
                    iv = find(Y(iYp:N) <= (Yp-hYpmin),1,'first') +iYp -1;
                    vlls(1,1) = iv;
                    vlls(1,2) = Y(iv);
                end
         end
         
            [Ymin0,iYmin0] =  min(Y(1:iYp(1)));
            iYv = [iYmin0; vlls(:,1)]; 
            Yv  = [ Ymin0; vlls(:,2)];
            V   = length(Yv);
            
            if P == V
                   [ymin iymin] = min (Y(iYp(P):N));
                   iYv = [iYv; iymin + iYp(P)- 1];
                   Yv  = [Yv; ymin];
            end
            
        if isempty(nYp)
            nYp = P;
        end

        if nYp < P
                  %sortYp  = 'descend';
                [~,js]    = sort(Yp,'descend');
                      jp  = sort(js(1:nYp),'ascend');  
                
                jv       = zeros(nYp+1,1);
                jv(1)    =  find(Yv == min(Yv(1:jp(1)  ,1)),1,'last');                
                if nYp>1
                    for    i = 2:nYp;
                               ind = iYp(jp(i-1))< iYv & iYv < iYp(jp(i));                       
                               jv(i) = find(Yv == min(Yv(ind,1)),1,'last');
                    end
                end
                jv(nYp+1)=  find(Yv == min(Yv(jp(nYp)+1:P+1,1)),1,'first');
                
                
                Yp   =  Yp(jp,1);
                iYp  = iYp(jp,1);
                
                Yv   = Yv (jv,1);
                iYv  = iYv(jv,1);
                
                P    = nYp;
                V    = P+1; 
        end
            
         % Determine Onset and Offset of every peak Yp (Matlab 2015 findpeaks)
         OnYp  = zeros(P,1);
         OffYp = zeros(P,1);
         
         for i = 1:P
             if i == 1
                   OnYp(i)  = iYv(1);
             else
                    jp = find(  Yp(1:i-1) >= Yp(i),1,'last');
                    if ~isempty(jp)
                       OnYp(i)  = iYv(jp+1);
                    else
                       OnYp(i)  = iYv(1);
                    end
             end
             
                if i== P
                    OffYp(i) = iYv(P+1);
                else                    
                    jp = find( Yp(i+1:P) >= Yp(i),1,'first');
                    if ~isempty(jp)
                       OffYp(i)  = iYv(jp+i);
                    else
                       OffYp(i)  = iYv(P+1);
                    end
                end
         end
         
               hYp = zeros(P,1);
              wYp  = zeros(P,1);
              dYp  = zeros(P,1);
         
if nargout > 2 
        dYp = [OnYp, OffYp];
        [hYp0,ihYp0] = max(Y(dYp),[],2);
        hYp =  Yp - hYp0;
end

if nargout >3 
        Onset   = zeros(P,1);
        Offset  = zeros(P,1);
          switch dYpType
                case 'height'                
                   iYl           = iYv(1:P);
                   iYr           = iYv(2:P+1);
                    
                    for i = 1:P
                        theta = Yp(i) - hYp(i)*wYpLevel;
                        iOn  =  find( Y(iYl(i):iYp(i))  <= theta,1,'last');
                        if isempty(iOn)
                            iOn = 1;
                        end
                        iOff =  find( Y(iYp(i):iYr(i))  <= theta,1,'first');
                        if isempty(iOff)
                            iOff = iYr(i)-iYp(i)+1;
                        end
                            Onset(i)  = iOn +iYl(i)-1;
                            Offset(i) =iOff +iYp(i)-1;
                    end
                   
                case 'prominence'
                    iYl           = ones(P,1)*1;
                    iYl(ihYp0==1) = dYp(ihYp0==1);
                    
                    iYr           = ones(P,1)*N;
                    iYr(ihYp0==2) = dYp(ihYp0==2);
                    for i = 1:P
                        theta = Yp(i) - hYp(i)*wYpLevel;
                            Onset(i)  = find( Y(iYl(i):iYp(i))  <= theta,1,'last') +iYl(i)-1;
                            Offset(i) = find( Y(iYp(i):iYr(i))  <= theta,1,'first')+iYp(i)-1;
                    end
          end          
       
        wYp = Offset-Onset;
end

if nargout > 4
        dYp = [Onset, Offset];
end

if nargout > 5
        y       = zeros(1,N);
        for i = 1:P
            y (OnYp(i):OffYp(i))= hYp0(i);
        end
            y(OffYp(P):N)       = hYp0(P);
        for i = 1:P
            y (Onset(i):Offset(i))= Yp(i);
        end
            
end
      
      switch sortYp
          case 'ascend'
              [~,jp]    = sort(Yp,'ascend');
              Yp  =  Yp(jp,:);
              iYp = iYp(jp,:);
              hYp = hYp(jp,:);
              wYp = wYp(jp,:);
              dYp = dYp(jp,:);
          case 'descend'
              [~,jp]    = sort(Yp,'descend');
              Yp  =  Yp(jp,:);
              iYp = iYp(jp,:);
              hYp = hYp(jp,:);
              wYp = wYp(jp,:);
              dYp = dYp(jp,:);
      end
end

