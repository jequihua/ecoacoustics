function S = TFRfilt(P,i0,j0,stft)
S = Ifilt(P,i0,j0,stft.amax,stft.amin,stft.thetha,stft.delta,stft.mode);

function I = Ifilt(P,i0,j0,Max,Min,alpha,delta,mode)
    I0 = zeros(size(P));
    P0 = P(i0,j0);
        area = 1;
        beta = 0;
    k = 0; 
    while area<Min | Max<area
    theta = db2pow(pow2db(P0)- alpha + beta);
    I1 = I0; 
    I1 (P>theta) = 1;
    if mode
        I = getconectset (I1,i0,j0,P);
    else
        I = getconectset (I1,i0,j0);
    end
    area = nnz(I)/nzmax(I);
        if Max<area
        beta = beta + delta;
        disp([beta area])
        elseif area<Min
        beta = beta - delta;
         disp([beta area])
        end

        k = k+1;
        if k ==50
            delta = delta/10;
        end
        if k >100 && Max<area
            break
        end    
    end
      function s = getconectset (S,i0,j0,varargin)
    jend = size(S,2);
    iend = size(S,1);
    if nargin>3
        X = varargin{1};
        mode = 2;
    else
        mode = 1;
    end
    
    for j = j0:jend
    i = find(S(:,j));
        if ~isempty(intersect(i,i0))  
            inull  = find(~S(:,j));
            iup = inull(find(inull>max(i0),1,'first'));
            idown = inull(find(inull<min(i0),1,'last'));
            if isempty(iup)
                iup = iend;
            end
            if isempty(idown)
                idown = 1;
            end
            S(iup:iend,j)= 0;
            S(1:idown,j) = 0;
            switch mode 
                case 1
                i0 = find(S(:,j));
                case 2
                    iX = find(S(:,j));
                    [v imax] = max(X(iX,j));
                    i0 = iX(imax);
            end
        else
            S(:,j) = 0;
            i0 = [];
        end
    end
    i0 = find(S(:,j0));
    for j = j0-1:-1:1
    i = find(S(:,j));
        if ~isempty(intersect(i,i0))  
            inull  = find(~S(:,j));
            iup = inull(find(inull>max(i0),1,'first'));
            idown = inull(find(inull<min(i0),1,'last'));
            if isempty(iup)
                iup = iend;
            end
            if isempty(idown)
                idown = 1;
            end
            S(iup:iend,j)= 0;
            S(1:idown,j) = 0;
            switch mode 
                case 1
                i0 = find(S(:,j));
                case 2
                   iX = find(S(:,j));
                    [v imax] = max(X(iX,j));
                    i0 = iX(imax);
            end
        else
            S(:,j) = 0;
            i0 = [];
        end
    end
    
    s = S;