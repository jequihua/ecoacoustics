function [e,f,t] = EFTcontour (E,F,T,A)
    
    
    J = size(E,2);
    
    EA = E; 
    EA(~A) = nan;
    
    e = nan(1,J);
    f = nan(1,J);
    t = nan(1,J);
    
    EF = EA.*(F*ones(1,size(EA,2)));
    f = nansum(EF)./nansum(EA);
        
   
    j = ~isnan(f);
    t(j) = T(j);   
    
   for j = 1:J; 
       if ~isnan(f(j))
        [v i] = min(abs(F-f(j)));        
        e(j) = E(i,j);        
       end
   end
   
   i = ~any(A,1);
   e(i)=[];
   f(i)=[];
   t(i)=[];