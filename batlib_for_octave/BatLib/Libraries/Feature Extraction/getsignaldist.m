function  [S,X] = getsignaldist(P,i0,j0,stft,mode)
if mode == 0
    S = Ifilt(P,i0,j0,stft.amax,stft.amin,stft.thetha,stft.delta,stft.mode);
     X =[];
elseif mode == 1
    [I,J] = size(P);
    A = I*J;
    iP0   = sub2ind([I J],i0,j0);
    
    P0  = pow2db(P(iP0));
    Beta  = stft.thetha;
    Delta = 0;
    theta = db2pow(P0-Beta+Delta);
    
    prop = stft.areatype;    
    X = Ifilt2(P,I,J,iP0,theta,prop,stft.conn);
    area = X.Area/A;
    
    if area>stft.amax
        Delta = stft.delta;
        Theta = (P0-Beta+Delta):Delta:P0;
        
        for theta = db2pow(Theta);
             X = Ifilt2(P,I,J,iP0,theta,prop,stft.conn);
             area = X.Area/A;
            if area<stft.amax
                break       
            end
        end
    end
    
    S = X.S;
    
    if nargout >1
        X.RefPoint = iP0;        
    end
    
end


    function X = Ifilt2(Z,I,J,iP0,theta,porperty,conn)  
    
    S0 = Z>theta;    
    X = struct();
    X.Connectivity  = conn;
    X.ImageSize     = [I,J];
    X.NumObjects    = 1;
    X.PixelIdxList  = cell(1,1);
    x = bwconncomp(S0,conn);    
        for k = 1:x.NumObjects
            if any(x.PixelIdxList{k} == iP0)
                    X.PixelIdxList{1}= x.PixelIdxList{k};                    
            end
        end
    area  = [porperty 'Area'];
    image = [porperty 'Image'];
    Y = regionprops(X,image,area,'BoundingBox','Area');
    
    i1 = uint16(Y.BoundingBox(2));
    i2 = uint16(Y.BoundingBox(4) +i1-1);
    j1 = uint16(Y.BoundingBox(1));
    j2 = uint16(Y.BoundingBox(3) +j1-1);    
    S = false(I,J);
    S(i1:i2,j1:j2) = Y.(image);
    
    X.PixelIdxList0 = X.PixelIdxList;
    X.PixelIdxList  = {find(S)};
    X.S             = S;
    X.Area0         = Y.Area;
    X.Area          = Y.(area);    
    X.BoundingBox   = Y.BoundingBox;
    y = regionprops(S,Z,'PixelValues');
    X.PixelValues   = y.PixelValues;
    
    
function I = Ifilt(P,i0,j0,Max,Min,alpha,delta,mode)
    I0 = zeros(size(P));
    P0 = P(i0,j0);
        area = 1;
        beta = 0;
    k = 0; 
    while area<Min || Max<area
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
        %disp([beta area])
        elseif area<Min
        beta = beta - delta;
         %disp([beta area])
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
                    [~,imax] = max(X(iX,j));
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
                    [~,imax] = max(X(iX,j));
                    i0 = iX(imax);
            end
        else
            S(:,j) = 0;
            i0 = [];
        end
    end
    
    s = S;