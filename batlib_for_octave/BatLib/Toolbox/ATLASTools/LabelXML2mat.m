function labels = LableXML2mat( filename )
%LABLEXML2MAT Summary of this function goes here
%   Detailed explanation goes here
xml = parseXML(filename);
labels = struct('startTime','endTime','timeStamp','value','comment','type','text','classentity');% = zeros(size(xml.Children,2),8);

for i=1:size(xml.Children,2)
    labels(i).startTime=str2double(xml.Children(i).Children(1).Children.Data);
    labels(i).endTime=str2double(xml.Children(i).Children(2).Children.Data);
    labels(i).timeStamp=str2double(xml.Children(i).Children(3).Children.Data);
    labels(i).value=str2double(xml.Children(i).Children(4).Children.Data);
    if size(xml.Children(i).Children(5).Children,1)~=0
        labels(i).comment=xml.Children(i).Children(5).Children.Data;
    end
    if size(xml.Children(i).Children(6).Children,1)~=0
        labels(i).type=xml.Children(i).Children(6).Children.Data;
    end
    if size(xml.Children(i).Children(7).Children,1)~=0
        labels(i).text=xml.Children(i).Children(7).Children.Data;
    end
    if size(xml.Children(i).Children(8).Children,1)~=0
        labels(i).classentity=xml.Children(i).Children(8).Children.Data;
    end
end

end

