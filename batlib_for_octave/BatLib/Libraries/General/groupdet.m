function [Ngroup,Index,Interval] = groupdet(T,m,gmin,gmax,dTmax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if nargin<5
    dTmax = 100e-3;
end
if nargin<4
    gmax = 2; %0 - 1 -inf
end
if nargin<3
    gmin = 1/2; %0 - 1 -inf
end
if nargin<2
    m = 1;
end

  nmax  = nzmax(T);
  Ngroup     = zeros(nmax,1);
  Index = zeros(nmax,2);
  Interval = zeros(nmax,2); 
  

    
if nmax>=2
    i = 1;
    g = 1;
    n = 1;
    dT = nan(1,m);
    dT(m) = T(i+1)-T(i);
    Index(g,1)   = i;
    Interval(g,1)= T(i);
    
    for i = 2:nmax        
            dTi  = T(i)-T(i-1); 
            if dTi <= dTmax                
                dTmean = nanmean(dT);
                e = dTi/dTmean;
                if gmin <= e && e <=  gmax
                    n = n+1;
                    j    = mod(n,m)+1;
                    dT(j)= T(i)-T(i-1);                                     
                else
                Ngroup(g)         = n;
                Index(g,2)   = i-1;
                Interval(g,2)= T(i-1);
                
                n = 1;
                g = g+1;                
                Index(g,1)   = i;
                Interval(g,1)= T(i);
                    if i<nmax
                    dT = nan(1,m);                
                    dT(m) = T(i+1)-T(i);
                    end
                end 
            else
                Ngroup(g)    = n;
                Index(g,2)   = i-1;
                Interval(g,2)= T(i-1);
                
                n = 1;
                g = g+1;                
                Index(g,1)   = i;
                Interval(g,1)= T(i);
                    if i<nmax
                    dT = nan(1,m);                
                    dT(m) = T(i+1)-T(i);
                    end
            end
    end
    Ngroup(g)    = n;
    Index(g,2)   = nmax;
    Interval(g,2)= T(nmax);
    
    
    Ngroup  = Ngroup(1:g);
    Index   = Index(1:g,1:2);
    Interval= Interval(1:g,1:2);  
    
else
    Ngroup       = 1;
    Index   = [1 1];
    Interval= [T(1) T(1)];
end

