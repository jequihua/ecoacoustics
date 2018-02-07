classdef AtlasLabelXML
    %ATLASLABELXML Summary of this class goes here
    %   Detailed explanation goes here
    
    % requires parseXML.m
    
    properties
    end
    
    methods(Static)
        function labels = loadXML( filename )
            %LABLEXML2MAT Summary of this function goes here
            %   Detailed explanation goes here
            xml = parseXML(filename);
            labels = struct('startTime',{},'endTime',{},'timeStamp',{},'value',{},'comment',{},'type',{},'text',{},'classentity',{});
            c = 0;
            for i=1:length(xml.Children)
                if strcmp(xml.Children(i).Name, 'label')    % perform on label tags only
                    c = c+1;
                    for j=1:length(xml.Children(i).Children)
                        switch xml.Children(i).Children(j).Name
                            case 'starttime'
                                labels(c).startTime=str2double(xml.Children(i).Children(j).Children.Data);
                            case 'endtime'
                                labels(c).endTime=str2double(xml.Children(i).Children(j).Children.Data);
                            case 'timestamp'
                                labels(c).timeStamp=str2double(xml.Children(i).Children(j).Children.Data);
                            case 'value'
                                labels(c).value=str2double(xml.Children(i).Children(j).Children.Data);
                            case 'comment'
                                if size(xml.Children(i).Children(j).Children,1)~=0
                                    labels(c).comment=xml.Children(i).Children(j).Children.Data;
                                end
                            case 'type'
                                if size(xml.Children(i).Children(j).Children,1)~=0
                                    labels(c).type=xml.Children(i).Children(j).Children.Data;
                                end
                            case 'text'
                                if size(xml.Children(i).Children(j).Children,1)~=0
                                    labels(c).text=xml.Children(i).Children(j).Children.Data;
                                end
                            case 'classentity'
                                if size(xml.Children(i).Children(16).Children,1)~=0
                                    labels(c).classentity=xml.Children(i).Children(16).Children.Data;
                                end
                        end
                    end
                end
            end
        end

        function  LabelTrack = saveXML( labels,name,externalChange,className,fileName )
            %MAT2LABELXML Summary of this function goes here
            %   Detailed explanation goes here

            LabelTrack = com.mathworks.xml.XMLUtils.createDocument('LabelTrack');

            root = LabelTrack.getDocumentElement;
            root.setAttribute('classname',className);
            root.setAttribute('externalchange',num2str(externalChange));
            root.setAttribute('name',name);

            for i=1:size(labels,2)
                label = LabelTrack.createElement('label');

                startTime=LabelTrack.createElement('starttime');
                startTime.appendChild(LabelTrack.createTextNode(num2str(labels(i).startTime)));
                label.appendChild(startTime);
                endTime=LabelTrack.createElement('endtime');
                endTime.appendChild(LabelTrack.createTextNode(num2str(labels(i).endTime)));
                label.appendChild(endTime);
                timestamp=LabelTrack.createElement('timestamp');
                timestamp.appendChild(LabelTrack.createTextNode(num2str(labels(i).timeStamp)));
                label.appendChild(timestamp);
                value=LabelTrack.createElement('value');
                value.appendChild(LabelTrack.createTextNode(num2str(labels(i).value)));
                label.appendChild(value);
                comment=LabelTrack.createElement('comment');
                comment.appendChild(LabelTrack.createTextNode(labels(i).comment));
                label.appendChild(comment);
                type=LabelTrack.createElement('type');
                type.appendChild(LabelTrack.createTextNode(labels(i).type));
                label.appendChild(type);
                text=LabelTrack.createElement('text');
                text.appendChild(LabelTrack.createTextNode(labels(i).text));
                label.appendChild(text);
                classentity=LabelTrack.createElement('classentity');
                classentity.appendChild(LabelTrack.createTextNode(labels(i).classentity));
                label.appendChild(classentity);

                root.appendChild(label);
            end
            xmlwrite(fileName,LabelTrack);
        end
    end
end

