function X = Lfeatures (featlist,P,I,E,F,T,C)

functlist = {   @CallDuration;...
                @PeakFrequency;...
                @StartFrequency;...
                @HalfTimeFrequency;...
                @EndFrequency;...
                @TopFrequency;...
                @BottomFrequency;...
                @FrequencyBandwidth;...
                @EnergyDistribution;...
                @Slope;...
                @FDerivative;...
                @DescriptiveFeatures;...                
                @HuMoments;...
                @ShapeFactors;...
                @ShapeFactors2;...
                @SplineDescriptiveFeatures};
if nargin == 0;
        functstr = cellfun(@(x) x(), functlist, 'UniformOutput', false);
        X = realign(functstr);        
elseif nargin==1;
    [ix, ifun] = getindex(featlist,functlist);    
    X = {ix,ifun};
    
elseif nargin==7    
    if iscellstr(featlist) || isempty(featlist)
        [ix, ifun] = getindex(featlist,functlist);        
    else 
        ix   = featlist{1};
        ifun = featlist{2};
    end
    
    deltaT = T(2)-T(1);
    dFdT = gradient(F,deltaT);     
    
   x = cellfun(@(x) x(P,I,E,F,T,dFdT,C), functlist(ifun), 'UniformOutput', false);
   x = cell2mat(x');
   X = x(ix);
end
function [ix, ifun] = getindex(featlist,functlist)
    if isempty(featlist)
        Nx = sum(cellfun(@(x) size(x(),1),functlist,'UniformOutput', true));
        Nfun =length(functlist);
        ifun = true(Nfun,1);
        ix   = true(Nx,1);
    elseif iscellstr(featlist)
        funcell = cellfun(@(x) x(), functlist, 'UniformOutput', false);
        Nfeat = length(featlist);
        Nfun  = length(funcell);
        A = zeros(Nfeat,Nfun);
        k = 0;
        for i = 1:Nfun
            for j = 1:Nfeat
                a = strcmpi(featlist{j},funcell{i});
                if any(a);
                    A(j,i) = find(a)+k;
                end
            end
                if any(A(:,i))
                k = k+size(funcell{i},1);
                end
        end
        ifun = any(A,1);
        ix = nonzeros(A)';
    end        
function Y = realign(X)
        if ~iscellstr (X)
            k = 1;
            for i = 1:length(X);
                if iscellstr (X(i))
                        Y(k,1) = X(i);
                        k = k+1;
                else
                    for j = 1:length(X{i})
                        Y(k,1) = X{i}(j);
                        k = k+1;
                    end
                end
            end
        else
            Y = X;
        end
    function x = CallDuration (~,~,~,~,T,~,~)
        if nargin == 0
            x = 'Call Duration';
        else
            x = T(end) - T(1);
        end
    function x = PeakFrequency (~,~,E,F,~,~,~)
        if nargin == 0
            x = 'Peak Frequency';
        else
            [~, iEmax] = max(E);
            x = F(iEmax);
        end
    function x = StartFrequency(~,~,~,F,~,~,~)
        if nargin == 0
            x = 'Start Frequency';
        else
            x = F(1);
        end
    function x = HalfTimeFrequency(~,~,~,F,~,~,~)
        if nargin == 0
            x = 'Frequency at half time';
        else
            x = F(ceil(length(F)/2));
        end            
    function x = EndFrequency(~,~,~,F,~,~,~)
        if nargin == 0
            x = 'End Frequency';
        else
            x   = F(end);
        end
    function x = TopFrequency (~,~,~,F,~,~,~)
        if nargin == 0
            x = 'Top Frequency';
        else
            x   = max(F);
        end
    function x = BottomFrequency (~,~,~,F,~,~,~)
        if nargin == 0
            x = 'Bottom Frequency';
        else
            x = min(F);
        end
    function x = FrequencyBandwidth (~,~,~,F,~,~,~)
        if nargin == 0
            x = 'Bandwidth';
        else
            x = max(F) - min(F);
        end
    function x = EnergyDistribution (~,~,E,~,~,~,~)
         Q = 4;     
        if nargin == 0
            x = cell(Q,1);
            strhead = 'Relative Energy per Quartile '; 
            for q = 1:Q;
            x{q} = [strhead num2str(q)];
            end
        else
            e = (E - min(E))./(max(E)- min(E));
            e = e./sum(e);
           
            q = quantile(e,0:1/Q:1);
            x = nan(1,Q);
            for i = 1:Q 
                ie = q(i)<= e & e<= q(i+1);
                x(i) = sum(e(ie));
            end
        end
    function x = Slope (~,~,~,F,T,~,~)
            S = 5;
           if nargin == 0
                x = cell(S,1);
                strhead = 'Slope ';                
                for s = 1:S;
                x{s} = [strhead num2str(s)];
                end
           else
                x = zeros(1,S);
                dT = (T(end)-T(1))/S;
                Tlim = T(1) + (0:S)*dT;
               for n = 1:S
                   iT = Tlim(n)<= T & T <=Tlim(n+1);
                   t = T(iT);
                   f = F(iT);
                   % Simple linear regression
                       alpha = sum((t-mean(t)).*(f-mean(f)));
                       beta  = sum((t-mean(t)).^2);           
                   x(n)= alpha/beta;
               end
           end
    function x = FDerivative(~,~,~,~,~,dFdT,~)
        if nargin == 0
            str = 'F-Derivative ';
            x = {   [str 'Absolute Sum'];...
                    [str 'Mean'];...
                    [str 'Absolute Mean'];...
                    [str 'Median'];...
                    [str 'Inflections']};
        else
            dFdTsign    = sign(dFdT);
            zrtest      = prod([dFdTsign(1:end-1); dFdTsign(2:end)])==-1;
            
            x = zeros(1,5);          
            x(1) = sum(abs(dFdT));
            x(2) = mean(dFdT);
            x(3) = mean(abs(dFdT));
            x(4) = median(dFdT);
            x(5)   = sum(zrtest);
        end
    function x = DescriptiveFeatures (~,~,~,F,T,dFdT,~)
            if nargin == 0                
            x = {   'Characteristic Frequency Time';...
                    'Characteristic Frequency';...
                    'Characteristic Frequency Slope';...
                    'Time of the Knee';...
                    'Frequency of the Knee';...
                    'Slope of the Knee';...
                    'Time of the 2nd Knee';...
                    'Frequency of the 2nd Knee';...
                    'Slope of the 2nd Knee'};
            else
                x = zeros(1,9);
                p = 0.6;
                dT = T(end)-T(1);
                iT = T <= (T(1) + p*dT);
                Df = abs(dFdT);
                iTc = find(Df == min(Df(iT)) & iT,1,'last');
                x(1) = T(iTc);
                x(2) = F(iTc);
                x(3) = dFdT(iTc);
                
                Df2 = abs(gradient(dFdT,dT));
                iT  = T <= T(iTc);
                iTk = find(Df2 == max(Df2(iT)) & iT,1,'first');                
                x(4) = T(iTk);
                x(5) = F(iTk);
                x(6) = dFdT(iTk);
                
                iT  = T > T(iTc);
                iTk = find(Df2 == max(Df2(iT)) & iT,1,'first');                
                x(7) = T(iTk);
                x(8) = F(iTk);
                x(9) = dFdT(iTk);
            end
    function x = CharFrequency(~,~,~,F,T,dFdT,~)
            if nargin == 0                
            x = {'Characteristic Frequency Time';...
                 'Characteristic Frequency';...
                 'Characteristic Frequency Slope'};
            else                
                p = 1 - 0.4;
                dT = T(end)-T(1);
                iT = T <= (T(1) + p*dT);
                Df = abs(dFdT);
                iTc = find(Df == min(Df(iT)) & iT,1,'first');
                x(1) = T(iTc);
                x(2) = F(iTc);
                x(3) = dFdT(iTc);
            end
    function x = KneeFrequency(~,~,~,F,T,dFdT,~)
            if nargin == 0
            x = {'Time of the knee';...
                 'Frequency of the knee'};
            else
                p = 1 - 0.4;                
                dT = T(end)-T(1);
                iT = T <= (T(1) + p*dT);                
                Df2 = abs(gradient(dFdT,dT));
                iTk = find(Df2 == max(Df2(iT)) & iT,1,'last');                
                x(1) = T(iTk);
                x(2) = F(iTk);               
            end
    function x = HuMoments(P,S,~,~,~,~,~)
        if nargin == 0
            str = ' Hu-Moment';
            x = {   ['1st' str];...
                    ['2nd' str];...
                    ['3rd' str];...
                    ['4th' str];...
                    ['5th' str];...
                    ['6th' str];...
                    ['7th' str]};
                    
        else
            p = P-min(P(:));
            I = S.*p;
            x = Hu_moments(I);            
        end
    function x = ShapeFactors(~,I,~,~,~,~,~)
        if nargin == 0
            str = ' Shape Factor';
            x = {   ['Aspect Ratio' str];...
                    ['Circularity' str];...
                    ['Enlongation' str];...
                    ['Eccentricity' str];...
                    ['Theta' str];...
                    ['Compactness' str]};
                    
        else
            Edge      = edge(I);         
            Area      = sum(sum(I));
            Perimeter = double(sum(sum(Edge)));
            
            m20 = mu_moment(2,0,I);
            m02 = mu_moment(0,2,I);
            m11 = mu_moment(1,1,I);
            
            A = 0.5 * (m20 + m02);
            B = 0.5 * sqrt(4 * m11^2 + (m20 - m02)^2);
            lamda1 = max([A - B, A + B]);       
            lamda2 = min([A - B, A + B]);
            [d_min, d_max] = max_min_diameter(I);
            
            AspectRatio = d_min / d_max;            
            Circularity = 4 * pi * Area / (Perimeter^2);
            Eccentricity= sqrt(1 - lamda2 / lamda1);
            Enlongation = sqrt(lamda2 / lamda1);
            Theta       = 0.5 * atan(2 * m11 / (m20 - m02));
            Compactness = (Area^2) / (2 * pi * sqrt(lamda2^2 + lamda1^2));
            
            x = [AspectRatio,...
                 Circularity,...
                 Enlongation,...
                 Eccentricity,...
                 Theta,...
                 Compactness];            
        end
    function x = ShapeFactors2(~,I,~,~,~,~,~)
        if nargin == 0
            str = ' Shape Factor';
            x = {   ['Orientation' str];...
                    ['Solidity' str];...
                    ['EquivDiameter' str];...
                    ['MajorAxisLength' str];...
                    ['MinorAxisLength' str];...
                    };
        else
            X = regionprops(I,'Orientation','Solidity','EquivDiameter','MajorAxisLength','MinorAxisLength');
            x = [X.Orientation,...
                 X.Solidity,...
                 X.EquivDiameter,...
                 X.MajorAxisLength,...
                 X.MinorAxisLength];            
        end
    function x = SplineDescriptiveFeatures(~,~,E,~,T,~,C)
                if nargin == 0                
                x = {   'Characteristic Frequency Time (Spline)';...
                        'Characteristic Frequency(Spline)';...
                        'Characteristic Frequency Slope (Spline)';...
                        'Time of the Knee (Spline)';...
                        'Frequency of the Knee (Spline)';...
                        'Slope of the Knee (Spline)';...
                        'Time of the 2nd Knee (Spline)';...
                        'Frequency of the 2nd Knee (Spline)';...
                        'Slope of the 2nd Knee (Spline)';...
                        'Number of Inflections (Spline)'};
                else
                    [~, iEmax] = max(E);
                    %Fp    = F(iEmax);
                    Tp    = T(iEmax);
                    
                    F0    = C.Fc;
                    T0    = C.Tc;
                    p     = C.PP.pieces; 
                    A     = C.PP.coefs;                    
                    k     = C.PP.order;
                    d     = C.PP.dim;
                    dA    = repmat(k-1:-1:1,d*p,1).* A(:,1:k-1);
                    ddA   = repmat(k-2:-1:1,d*p,1).*dA(:,1:k-2);                    
                    dPP   = mkpp(T0,dA,d);
                    ddPP  = mkpp(T0,ddA,d);
                    
                    dF0   = ppval(dPP,T0);
                    ddF0  = ppval(ddPP,T0);
                    zdF   = sign(dF0(1:p)).*sign(dF0(2:p+1))<=0;
                    zddF   = sign(ddF0(1:p)).*sign(ddF0(2:p+1))<=0;
                    if any(zdF);
                        dF0int = nan(1,p+1); 
                        for z = find(zdF);
                            if zddF(z)
                                tm       = -ddA(z,2)/ddA(z,1)+T0(z);                                
                                mdF(1)   = dA(z,1)*(tm^2) + dA(z,2)*tm +dA(z,3);
                                mdF(2:3) =  dF0(z:z+1);                                
                                dF0int(z) = abs(max(mdF) - min(mdF));
                            else
                                dF0int(z) = abs(dF0(z) - dF0(z+1));
                            end
                        end
                        [~,ic] = min(dF0int);                            
                        Tr  = roots(dA(ic,:))+T0(ic);                        
                        Tc  = Tr(T0(ic)<= Tr & Tr <T0(ic+1));
                        Fc  = polyval(A(ic,:),Tc-T0(ic));
                        dFc = polyval(dA(ic,:),Tc-T0(ic));                        
                    else
                        [dFc,ic] = min(abs(dF0)); 
                        Fc = F0(ic);
                        Tc = T0(ic);                     
                    end
                    
                        [~,ik1] = max(abs(ddF0(1:ic)));
                        Fk1 = F0(ik1);
                        Tk1 = T0(ik1);
                        dFk1 = dF0(ik1);
                        
                        [~,ik2] = max(abs(ddF0(ic:end)));
                        ik2 = ik2+ic-1;
                        Fk2 = F0(ik2);
                        Tk2 = T0(ik2);
                        dFk2 = dF0(ik2);
                    
                    x = [Tc-Tp,Fc,dFc,...
                         Tk1-Tp,Fk1,dFk1,...
                         Tk2-Tp,Fk2,dFk2,...
                         sum(zdF)];
                    
                end