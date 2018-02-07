function [e,f,t,C] = EFTcontour2 (E,F,T,S,delta)
    
    E(~S) = nan;
    i = any(E,2);
    j = any(E,1);
    I = nnz(i);
    J = nnz(j);
    
if nargin ==4  
    EF = E(i,j).*repmat(F(i),1,J);
    f = nansum(EF,1)./nansum(E(i,j),1);    
    t = T(j);
    e = interpE(E(i,j),F(i),f);
    C    = [];
elseif nargin == 5
 
    dT   = T(2)-T(1);
    minT = min(T(j));
    maxT = max(T(j));
    durT = maxT-minT;
    if delta < 1        
        Tdelta = max(2*dT,delta);
    elseif delta >= 1
        Tdelta = max(2*dT,durT/delta);
    end
   

    K = ceil(durT/Tdelta);
   
    if K < 4
        delta = 4;
        Tdelta = max(2*dT,durT/delta);
        K = ceil(durT/Tdelta);
    end
    
    ec = nan(1,K);
    tc = nan(1,K);
    fc = nan(1,K);
    ts = nan(1,K);
    fs = nan(1,K);
    
    for k  = 1:K        
    jk   = minT +(k-1)*Tdelta <= T & T <= minT + k*Tdelta;
    Jk   = nnz(jk); 
    Ek   = E(i,jk);
    ek   = nansum(Ek(:));
    
    FF = repmat(F(i),1,Jk);
    TT = repmat(T(jk),I,1);
    
    ec(k) = nanmean(Ek(:)); 
    fc(k) = sum(nansum(Ek.*FF))/ek;
    tc(k) = sum(nansum(Ek.*TT))/ek;
    fs(k) = sum(nansum( ((FF-fc(k)).^2).*Ek))/ek;
    ts(k) = sum(nansum( ((TT-tc(k)).^2).*Ek))/ek;
    end
    
    pp = spline(tc,fc);
    
    t = T(j);
    f = ppval(pp,t);
    e = spline(tc,ec,t);
    
    C.Tc  = tc;
    C.Fc  = fc;
    C.Tvar= ts;
    C.Fvar = fs;
    C.PP   = pp;
end

function e = interpE(E,F,f)
    e = nan(size(f));
    for j = 1:length(f);
        [~,iE] = min(abs(F-f(j)));        
        e(j) = E(iE,j);        
    end
