function [iA,B] = interval (A,I,varargin)
if length(I)~=2
    error('myApp:argChk', 'Wrong number of input arguments')
end

if isempty(varargin)
    type = '[]';
else
    type = varargin{1};
end

switch type
    case '[]'
    iA = I(1)<= A & A<=I(2);
    case '()'
    iA = I(1)<  A & A< I(2);
    case '(]'
    iA = I(1)<  A & A<=I(2);
    case '[)'
    iA = I(1)<= A & A< I(2);
end

if nargout == 2
B = A(iA);
end
        