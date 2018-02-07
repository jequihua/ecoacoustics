function [D,Ex,Dx,Px,iPx,Ey,Dy,Py,iPy] = PdetectorYX(Z,alpha,beta,theta)
%Detection per maximum energy bands

    [I,J] = size(Z);
        D = false(I,J);
    
    e = pow2db(nansum(Z,2));    
    if theta (1)
        Ey = e;
    else
        Ey  = smooth(e,theta(1));
    end
    
    if (max(Ey)- min(Ey))>= alpha(1)
            
        [Dy, idy1, idy2,Py,ipy] = Pdetector1(Ey,alpha(1),beta(1));  
        iPy = [idy1 ipy idy2];
        Ky   = length(ipy);
                
        Dx  = false(Ky,J);
        Px  = cell(Ky,1);
        iPx = cell(Ky,1);        
        Ex   = nan(Ky,J);
    
        for ky = 1:Ky
            i = iPy(ky,1):iPy(ky,3);
            z = Z(i,:);
            e = pow2db(nansum(z,1));            
            if theta(1)               
                Ex (ky,:) = e;
            else
                Ex (ky,:) = smooth(e,theta(2));
            end            
        
            e = Ex(ky,:);            
            if (max(e)- min(e))>= alpha(2)
                [dx, idx1, idx2,px,ipx] = Pdetector1(e,alpha(2),beta(2));                
                if ipx
                    Dx(ky,:)   = dx;
                    Px{ky,1}   = px;
                    iPx{ky,1} = [idx1 ipx idx2];                    
                    
                    Kx    = nnz(ipx);
                    for kx = 1:Kx
                        D(ipy(ky),ipx(kx)) = 1;
                    end                
                end
            end
        end
        
    end
        
end

