function X = Gfeatures (featlist,Data)
    functlist = {   @InterPulseDuration;...
                    @PeakFrequencyAlternation;...
                    @FDerivateAlternation};
    if nargin == 0;
        functstr = cellfun(@(x) x(), functlist, 'UniformOutput', false);
        X = realign(functstr);        
    elseif nargin==1;
        
    [ix ifun] = getindex(featlist,functlist);    
    X = {ix,ifun};
    
    elseif nargin==2    
    if iscellstr(featlist) || isempty(featlist)
        [ix ifun] = getindex(featlist,functlist);        
    else 
        ix   = featlist{1};
        ifun = featlist{2};
    end   
   x = cellfun(@(x) x(Data), functlist(ifun), 'UniformOutput', false);
   x = cell2mat(x');
   X = x(:,ix);
    end
function [ix ifun] = getindex(featlist,functlist)
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
    function x = InterPulseDuration(Data)
        if nargin == 0
            x = {'Interpulse Duration'};
        else
            x = nan(Data.N,1);
            for n = 1:Data.N-1
                t1 = Data.EFTcurve{n+1,3}(1) + Data.TimeLim2(n+1,1);
                t0 = Data.EFTcurve{n,3}(1)   + Data.TimeLim2(n,1);
                x(n) = t1-t0;
            end
                x(Data.N) = x(Data.N-1);
        end
    function x = PeakFrequencyAlternation(Data)
        if nargin == 0
            x = {'Peak Frequency Alternation'};
        else
            i = find(strcmp(Data.FeatListStr, 'Peak Frequency'));
            x = nan(Data.N,1);
            for n = 1:Data.N-1
                x(n) = Data.Data(n+1,i) - Data.Data(n,i);
            end
                x(Data.N) = x(Data.N-1);
        end
        
    function x = FDerivateAlternation(Data)
        if nargin == 0
            x = {   'F-Derivative Absolute Sum Alternation';...
                    'F-Derivative Mean Alternation';...
                    'F-Derivative Absolute Mean Alternation';...
                    'F-Derivative Median Alternation';...
                    'F-Derivative Inflections Alternation'};
        else
            featlist = {'F-Derivative Absolute Sum';...
                        'F-Derivative Mean';...
                        'F-Derivative Absolute Mean';...
                        'F-Derivative Median';...
                        'F-Derivative Inflections'};
            I = size(featlist,1);
            k = zeros(1,I);
            for i = 1:I
                k(i) = find(strcmp(Data.FeatListStr, featlist{i}));
            end
            
            x = nan(Data.N,1);
            for i = 1:I
                for n = 1:Data.N-1
                    x(n,i) = Data.Data(n+1,k(i)) - Data.Data(n,k(i));
                end
                    x(Data.N) = x(Data.N-1);
            end
        end

