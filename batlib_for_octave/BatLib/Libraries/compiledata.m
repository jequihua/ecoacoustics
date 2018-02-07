function compiledata(data)

indexfield = {'FileIndex','DetIndex'};
audiofield = {'Name','Path','TimeLength'};
 detdfield   = {'DetectionTimePoint','DetectionFreqPoint','DetectionTimeLength','NumberOfHarmonics'};%,...
%                 'FrequencyContourMean','FrequencyContourStd',...
%                 'FrequencyMaxMean','FrequencyMaxStd'};
%featfield   = {''};
N  = data.audioinfo.FilesAmount;
X0 = cell(N,9);
X0(1,:) = {'FileIndex','Path','Name','File Duration'...
            'NumberOfDetections',...
            'DetectionTimeLengthMean','DetectionTimeLengthStd',...    
            'FrequencyPeakMean','FrequencyPeakStd'};
i = 2;
for n = find(data.audioinfo.Key)'
    X0{i,1} = data.audioinfo.FileIndex(n);
    X0{i,2} = data.audioinfo.Path{n};
    X0{i,3} = data.audioinfo.Name{n};
    X0{i,4} = data.audioinfo.TimeLength(n);
    X0{i,5} = data.detdata.NumberOfDetections(n);
    X0{i,6} = mean(data.detdata.DetectionTimeLength{n});
    X0{i,7} = std(data.detdata.DetectionTimeLength{n});
    X0{i,8} = mean(data.detdata.DetectionFreqPoint{n});
    X0{i,9} = std(data.detdata.DetectionFreqPoint{n});
    i=i+1;
end

    system    = data.path.Sytem;
    winheader = data.path.WinHeader;
    linheader = data.path.LinHeader;
    filepath  = data.path.Worksheets;
    filepath = fullfile( changepath(filepath,system,winheader,linheader),data.Title);
    filepath = [filepath '.xlsx'];

p = 1;
xlswrite(filepath,X0,p);

B = length(indexfield);
A = length(audiofield);
D = length(detdfield);
F=0;
% F = size(featfield,1);
I = sum(data.detdata.Key)+1;
J = B+A+D+F;

X = cell(I,J);
j = 1;

    for b = 1:B  
            X{1,j} = indexfield{b};
            j=j+1;
    end
    for a = 1:A  
            X{1,j} = audiofield{a};
            j=j+1;
    end
    for d = 1:D  
            X{1,j} = detdfield{d};
            j=j+1;
    end
%     for f = 1:F  
%             X{1,j} = featfield{f,2};
%             j=j+1;
%     end


i = 2;
k = 1;
N = nnz(data.audioinfo.Key);
for n = find(data.audioinfo.Key)'
    s = 1;
    for m = 1:data.detdata.NumberOfDetections(n)
        j = 1;
        for b = 1:B
            x = {num2str(k);num2str(s)};
            X(i,j) = x(b);
            j=j+1;
        end
        for a = 1:A
                x = data.audioinfo.(audiofield{a})(n);
            if iscellstr(x) 
                X(i,j) = x;
            else
                X{i,j} = x;
            end
            j=j+1;
        end
        for d = 1:D
                x = data.detdata.(detdfield{d}){n}(m);
            if iscellstr(x) 
                X(i,j) = x;
            else
                X{i,j} = x;
            end
            j=j+1;
        end
%         for f = 1:F
%                 x = data.featext.(featfield{f,1}).(featfield{f,2}){n}(s);
%             if iscellstr(x) 
%                 X(i,j) = x;
%             else
%                 X{i,j} = x;
%             end
%             j=j+1;
%         end
        s = s+1;
        i = i+1;
    end
        k = k+1;
        disp(['Compiling file: ' num2str(n) '/' num2str(N) ': ' data.audioinfo.Name{n}])
end



maxI = 65500;
if I> maxI;
        i = [1:maxI:I,I];
    for k = 1:length(i)-1
    Xk = [X(1,:); X( i(k):i(k+1) , : )];    
    xlswrite(filepath,Xk,p+k)
    end;
else
    xlswrite(filepath,X,p+1)
end

