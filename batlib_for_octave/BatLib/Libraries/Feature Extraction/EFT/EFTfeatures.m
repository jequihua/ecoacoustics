function X = EFTfeatures (featlist,E,F,T)

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
                @CharFrequency;...
                @KneeFrequency};
if nargin == 0;
        functstr = cellfun(@(x) x(), functlist, 'UniformOutput', false);
        X = realign(functstr);        
elseif nargin==1;
    [ix, ifun] = getindex(featlist,functlist);    
    X = {ix,ifun};
    
elseif nargin==4    
    if iscellstr(featlist) || isempty(featlist)
        [ix, ifun] = getindex(featlist,functlist);        
    else 
        ix   = featlist{1};
        ifun = featlist{2};
    end
    
    deltaT = T(2)-T(1);
    dFdT = gradient(F,deltaT);     
    
   x = cellfun(@(x) x(E,F,T,dFdT), functlist(ifun), 'UniformOutput', false);
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
    function x = CallDuration (~,~,T,~)
        if nargin == 0
            x = 'Call Duration';
        else
            x = T(end) - T(1);
        end
    function x = PeakFrequency (E,F,~,~)
        if nargin == 0
            x = 'Peak Frequency';
        else
            [~, iEmax] = max(E);
            x = F(iEmax);
        end
    function x = StartFrequency(~,F,~,~)
        if nargin == 0
            x = 'Start Frequency';
        else
            x = F(1);
        end
    function x = HalfTimeFrequency(~,F,~,~)
        if nargin == 0
            x = 'Frequency at half time';
        else
            x = F(ceil(length(F)/2));
        end            
    function x = EndFrequency(~,F,~,~)
        if nargin == 0
            x = 'End Frequency';
        else
            x   = F(end);
        end
    function x = TopFrequency (~,F,~,~)
        if nargin == 0
            x = 'Top Frequency';
        else
            x   = max(F);
        end
    function x = BottomFrequency (~,F,~,~)
        if nargin == 0
            x = 'Bottom Frequency';
        else
            x = min(F);
        end
    function x = FrequencyBandwidth (~,F,~,~)
        if nargin == 0
            x = 'Bandwidth';
        else
            x = max(F) - min(F);
        end
    function x = EnergyDistribution (E,~,~,~)
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
    function x = Slope (~,F,T,~)
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
    function x = FDerivative(~,~,~,dFdT)
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
        function x = CharFrequency(~,F,T,dFdT)
            if nargin == 0                
            x = {'Characteristic Frequency Time';...
                 'Characteristic Frequency';...
                 'Characteristic Frequency Slope'};
            else                
                p = 1- 0.4;
                dT = T(end)-T(1);
                iT = T >= T(1) + p*dT;
                Df = abs(dFdT);
                iTc = find(Df == min(Df(iT)) & iT,1,'first');
                x(1) = T(iTc);
                x(2) = F(iTc);
                x(3) = dFdT(iTc);
            end
        function x = KneeFrequency(~,F,T,dFdT)
            if nargin == 0
            x = {'Time of the knee';...
                 'Frequency of the knee'};
            else
                p = 1- 0.4;                
                dT = T(2)-T(1);
                iT = T < T(1) + p*dT;                
                Df2 = abs(gradient(dFdT,dT));
                iTk = find(Df2 == max(Df2(iT)) & iT,1,'last');                
                x(1) = T(iTk);
                x(2) = F(iTk);               
            end