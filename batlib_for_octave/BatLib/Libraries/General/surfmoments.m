function [Mx,My,C,I] = surfmoments(x,y,Z)


[X,Y] = meshgrid(x,y);
 
Ex  = sum(Z,1);
Ey  = sum(Z,2);
E   = sum(Ex);

ZX  = Z.*X;
ZY  = Z.*Y;

Cx = sum(ZX(:))/E;
Cy = sum(ZY(:))/E;

mx = sum(ZX,1)./Ex;
my = sum(ZY,2)./Ey;
[MX,MY]  = meshgrid (mx,my);

sx = sum( ((X-MX).^2).*Z , 1)./Ex; 
sy = sum( ((Y-MY).^2).*Z , 2)./Ey;

Ix = sum(sum( ((X-Cx).^2).*Z))/E;
Iy = sum(sum( ((Y-Cy).^2).*Z))/E;

Mx = [mx; sx];
My = [my, sy];

C  = [Cx,Cy];
I  = [Ix,Iy];
    
 