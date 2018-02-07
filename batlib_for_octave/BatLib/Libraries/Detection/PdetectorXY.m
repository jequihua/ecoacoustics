function [D,Ex,Dx,Px,iPx,Ey,Dy,Py,iPy] = PdetectorXY(Z,alpha,beta,delta,scale)
%Detection per maximum energy bands

    [I,J] = size(Z);
        D = false(I,J);
        
        
        
    if nargin == 4
        scale = 'dB';
    end
    
    switch scale
        case 'dB'         
        case 'pow'
         alpha = pow2db(-alpha);         
        case 'mag'
         alpha = mag2db(-alpha);
    end
    
        
        Dx  = false(1,J);
        Ex = sumZ(Z,1,delta(1),scale);
        Px  = [];
        iPx = [];       
    
    if (max(Ex)- min(Ex))>= alpha(1)            
        [Dx, idx1, idx2,Px,ipx] = Pdetector1(Ex,alpha(1),beta(1));  
        iPx = [idx1 ipx idx2];
        Kx   = length(ipx);
                
        Dy  = false(I,Kx);
        Ey   = nan(I,Kx);
        Py  = cell(1,Kx);
        iPy = cell(1,Kx);
    
        for kx = 1:Kx
            j = iPx(kx,1):iPx(kx,3);
            z = Z(:,j);
            Ey (:,kx) = sumZ(z,2,delta(2),scale);
            e = Ey(:,kx);            
            if (max(e)- min(e))>= alpha(2)
                [dy, idy1, idy2,py,ipy] = Pdetector1(e,alpha(2),beta(2));                
                if ipy
                    Dy(:,kx)   = dy;
                    Py{1,kx}   = py;
                    iPy{1,kx} = [idy1 ipy idy2];                    
                    
                    Ky    = nnz(ipy);
                    for ky = 1:Ky
                        D(ipy(ky),ipx(kx)) = 1;
                    end                
                end
            end
        end        
    else
        Dy  = [];
        Ey  = sumZ(Z,2,delta(2),scale);
        Py  = {};
        iPy = {};    
    end 
end

function e = sumZ(z,dim,delta,scale)
            switch scale
                case 'dB'
                 e     = pow2db( nansum(z,dim));
                case {'pow','mag'}
                 e     = nansum(z,dim);
            end
            
            if delta               
                e = smooth(e,delta);
            end
end
