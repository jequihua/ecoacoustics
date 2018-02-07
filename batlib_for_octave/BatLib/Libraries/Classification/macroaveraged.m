function [F,P,R] = macroaveraged(M,C)
    if nargin <2
    C = size(M,1);
    end
    p = diag(M)./sum(M,1)';
    r = diag(M)./sum(M,2);
    f = 2*(p.*r)./(p+r);
    
    P = nansum(p)/C;
    R = nansum(r)/C;
    F = nansum(f)/C;
    
    if isnan(P)
        P = 0;
    end
    if isnan(R)
        R = 0;
    end
    if isnan(F)
        F = 0;
    end