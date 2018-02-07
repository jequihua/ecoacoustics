function C = rgb2cnum(rgb,A)

    if nargin == 1 || isempty(A)
        Ax = 'FF';
    else
        Ax = num2str(A,'%02X');
    end
    Rx = num2str(rgb(1),'%02X');
    Gx = num2str(rgb(2),'%02X');
    Bx = num2str(rgb(3),'%02X');
    
    Cx = [Ax Rx Gx Bx];
    C = hex2dec(Cx)- 2^32;