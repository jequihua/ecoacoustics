function  haxes  = plotspec(x,y,Z,aspect,dBmode,h)
    
    if nargin == 5
       h = 1;
    end
    [i,j] = size(Z);
    if isempty(x)        
        x = 1:j;
    end
    if isempty(y)
        y = 1:i;
    end
        
        L = aspect(1);
        H = aspect(2);
        htable  = reshape(1:H*L,L,H)';
        hytable = htable(1:H-1,1);
        hxtable = htable(H,2:L);
        hztable = htable(1:H-1,2:L);
        
        if dBmode
            Ey = pow2db(nansum(Z,2));
            Ex = pow2db(nansum(Z,1));
            E  = pow2db(Z);
        else
            Ey = nansum(Z,2);
            Ex = nansum(Z,1);
            E = Z;
        end
        
        YLim = [min(y) max(y)];
        XLim = [min(x) max(x)];
        iE   = ~isinf(E); 
        minE = min(min(E(iE),[],1),[],2);
        maxE = max(max(E(iE),[],1),[],2);
        meanE = nanmean(nanmean(E(iE),1),2);
        
        Zlim = [minE max(maxE,1)];
        CLim = [meanE maxE];
        
        hy  = subplot (H,L,hytable);
        plot(hy,Ey,y,'-b');
        ylim(YLim)    
        
        
        hx  = subplot (H,L,hxtable);        
        plot(hx,x,Ex)
        xlim(XLim);
        
        
        hz  = subplot (H,L,hztable(:));        
        surf(hz,x,y,E,'edgecolor','none')
        axis tight, view(0,90)
        ylim(YLim);
        xlim(XLim);
        zlim(Zlim);
        caxis(CLim);
        
        haxes = [h hx hy hz];
end

