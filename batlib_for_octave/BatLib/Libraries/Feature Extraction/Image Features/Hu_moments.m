function M = Hu_moments(A)

% This function Calculates the Seven Invariant Moments for the image A
% the output of this function is a Vector M ; called the Feature vector
% the vector M is a column vector containing M1,M2,....M7

n20=nu_moment(2,0,A);
n02=nu_moment(0,2,A);
n11=nu_moment(1,1,A);
n30=nu_moment(3,0,A);
n12=nu_moment(1,2,A);
n21=nu_moment(2,1,A);
n03=nu_moment(0,3,A);

% First Moment
M1=n20+n02;

% Second Moment
M2=(n20-n02)^2 + 4*n11^2;

% Third Moment
M3=(n30-3*n12)^2 + (3*n21-n03)^2;

% Fourth Moment
M4=(n30+n12)^2 + (n21+n03)^2;

% Fifth Moment
M5 =(n30-3*n21)*(n30+n12)*((n30+n12)^2 - 3*(n21+n03)^2) +...
    (3*n21-n03)*(n21+n03)*(3*(n30+n12)^2 - (n21+n03)^2);

% Sixth Moment
M6=(n20-n02)*((n30+n12)^2-(n21+n03)^2) + 4*n11*(n30+n12)*(n21+n03);

% Seventh Moment
M7=(3*n21-n03)*(n30+n12)*((n30+n12)^2-3*(n21+n03)^2) - ...
   (n30+3*n12)*(n21+n03)*(3*(n30+n12)^2-(n21+n03)^2);

M=[M1    M2     M3    M4     M5    M6    M7];