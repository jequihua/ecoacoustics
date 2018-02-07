function [I,C] = isclipped(y,t,w,fs,Cmax)
I = false(size(t));
C = zeros(size(t));
s = t*fs+1;
smax = nzmax(y);


if nargin<5
    Cmax=3;
end

for i = 1:length(s)
    s1  = uint64(max(s(i)- w/2,1));
    s2  = uint64(min(s(i)+ w/2,smax));
    x = y(s1:s2);
    
    
    xmax = max(x);
    xmin = min(x);    
    
    N          = length(x);
    ixx        = false(1,N);
    ixx(1:N-1) = x(1:N-1)== x(2:N);
    
    Nxmax = sum(x(ixx)==xmax);
    Nxmin = sum(x(ixx)==xmin);   
    
    if Nxmax+Nxmin > Cmax
        I(i) = true;
        C(i) = (Nxmax+Nxmin)/N;
    end
    
end
    
    
    