function [X0,X1] = compiledata2(data,Xfile,Xcall)


I  = data.audioinfo.FilesAmount;
J  = size(Xfile,1);

X0 = cell(I+1,J);

for j = 1:J;
    X0{1,j} = Xfile{j,1};    
end

for i = 1:I;
    for j = 1:J
        try
        x = Xfile{j,2};
        f = Xfile{j,3};
        k = data.detdata.Selection.Key{i};
        X0{i+1,j} = f(x,i,k);
        end
    end
end


I = sum(data.detdata.Selection.Amount);
J = size(Xcall,1);
X1 = cell(I+1,J);

for j = 1:J;
    X1{1,j} = Xcall{j,1};    
end


i = 1;
for n = find(data.detdata.Selection.Amount)';
    k = find(data.detdata.Selection.Key{n})';
    M  = data.detdata.Selection.Amount(n);
    for m = 1:M;
        for j = 1:J
            x = Xcall{j,2};
            f = Xcall{j,3};
            if ~isempty(x)
                try
                X1{i+1,j} = f(x,n,k(m),m); 
                end
            end
        end
        i = i+1;
    end
end


system    = data.path.Sytem;
%     winheader = data.path.WinHeader;
%     linheader = data.path.LinHeader;
%     filepath  = data.path.Worksheets;
    filepath = fullfile( changepath(data.path.Worksheets,data.path),data.Title);
    %filepath = [filepath '.xlsx'];

p = 1;
    try
    sheetname0 = ' TABLE A (Audiofile Summary)';
    xlswrite([filepath '.xlsx'],X0,sheetname0);
    
    maxI = 65500;
    if I> maxI;
            i = [1:maxI:I,I];
        for k = 1:length(i)-1
        Xk = [X1(1,:); X1( i(k):i(k+1) , : )];    
        xlswrite([filepath '.xlsx'],Xk,sheetname1(k));        
        end;
    else
        xlswrite([filepath '.xlsx'],X1,sheetname1);    
    end
    end

    exptable(X0,[filepath sheetname0],'.txt')
    exptable(X1,[filepath sheetname1()],'.txt');
    
end


    function str  = sheetname1(page)
            if nargin == 0
             str = ' TABLE B (Calls Summary)';
            else
             str = ['TABLE B' num2str(page,'%2.0f') ' (Calls Summary)'];
            end
    end
    function exptable(X,filepath,ext)
    varnames = getvalidnames(X);
    T = array2table(X(2:end,:),'VariableNames',varnames);
    writetable(T,[filepath ext],'Delimiter','|');
    end
    function varnames = getvalidnames(X)
        I = size(X,2);
        varnames  = cell(1,I);
        for i = 1:I
            str   = X{1,i};
            if ~isvarname(str)
                if isempty(str)
                str = ['Var' num2str(i)];
                else
                    str   = str(~isspace(str));  
                    str   = str(isstrprop(str,'alphanum'));
                end
            end
            varnames{1,i} = str;
        end
    end
