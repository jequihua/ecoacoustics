function [result] = csv2toolMat(file)
fid = fopen(file);
data = textscan(fid,'%f,%f,%f');
fclose(fid);
resultCount = 1;
for i=1:size(data{1},1)
    result(1,resultCount)=data{1}(i);
    result(2,resultCount)=data{3}(i);
    resultCount = resultCount+1;
    for ii=100:100:floor((data{2}(i)-data{1}(i)))-100
        result(1,resultCount)=data{1}(i)+ii;
        result(2,resultCount)=data{3}(i);
        resultCount = resultCount+1;
    end
end

end
