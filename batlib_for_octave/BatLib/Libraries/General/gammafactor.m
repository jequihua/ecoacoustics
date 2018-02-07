function [gamma, dt] = gammafactor(t,m)

if nargin<1
    m =1;
end

    N     = length(t);
    gamma = ones(1,N);
    dt    = nan(1,N);
    
if N>=2
    dt(2:N) = diff(t);
    dt(1)   = dt(2);
    if N>=3          
        for i = 3:N
            gamma(i) = dt(i)/mean(dt(max(i-m,1):i-1));
        end
    end
else
    dt = [];
    gamma = [];
end