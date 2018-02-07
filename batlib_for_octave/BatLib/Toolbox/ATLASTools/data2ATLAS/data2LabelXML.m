function XMLlabel = data2LabelXML (data,nfile,ltrack,ProjectPath)
    
    Nltrack         = size(ltrack,2);
for i =  Nltrack
    
    filepath = fullfile(ProjectPath.labeltracks,ltrack{i}.Text);
    templatepath = fullfile(data.Atlas.TemplatePath,'labeltracks',ltrack{i}.Text);
    XMLlabel = xmlread(templatepath);
    
    switch ltrack{1}.Text
        case 'SpID.xml'
        N       = data.detdata.Selection.Amount(nfile);
        Index = data.detdata.Selection.DetIndex{nfile};
        DetTimeInt = data.detdata.VADmethod.DetectionTimeInt{nfile}(Index,:)*1e3;
        DetTimeDur = ceil(DetTimeInt(:,2)) - floor(DetTimeInt(:,1));
        t = DetTimeDur * [0 0.5 1];
        y = ones(size(t))*0.5;
        
        starttime = cell(1,N) ;
        endtime   = cell(1,N); 
        timestamp = cell(1,N);
        for n = 1:N
            starttime{n}     = num2str(floor(DetTimeInt(n,1)),'%1.0f');
            endtime{n}       = num2str(ceil (DetTimeInt(n,2)),'%1.0f');
            timestamp{n}     = num2str(round(now));
        end
        
        
        value(1:N)    = {'1.0'};
        comment(1:N)  = {''};
        Type(1:N)     = {'MANUAL'};
        Text          = cellstr(num2str(Index','Call%2.0f'));
        
        if isfield(data,'classinfo')                
        else
                classentity    =  cell(1,N);
            for n = 1:N
                classentity{n} = 'undefined';
            end
        end
    end
    
   
   
        name            = char(XMLlabel.getDocumentElement.getAttribute('name'));
        externalchange  = char(XMLlabel.getDocumentElement.getAttribute('externalchange'));
        classname       = char(XMLlabel.getDocumentElement.getAttribute('classname'));
        isContinuous    = char(XMLlabel.getDocumentElement.getAttribute('isContinuous'));
        interpolationType = char(XMLlabel.getDocumentElement.getAttribute('interpolationType'));
        
        LabelTrack = com.mathworks.xml.XMLUtils.createDocument('LabelTrack');    
            root = LabelTrack.getDocumentElement;
            root.setAttribute('name',name);
            root.setAttribute('classname',classname);
            root.setAttribute('externalchange',num2str(externalchange));
            root.setAttribute('isContinuous',isContinuous);
            root.setAttribute('interpolationType',interpolationType);
        
    
    for n = 1:N
        label = LabelTrack.createElement('label');       
        label  = appendChild(LabelTrack,label,'starttime',starttime{n});        
        label  = appendChild(LabelTrack,label,'endtime',endtime{n});
        label  = appendChild(LabelTrack,label,'timestamp',timestamp{n});
        label  = appendChild(LabelTrack,label,'value',value{n});
        label  = appendChild(LabelTrack,label,'comment',comment{n});
        label  = appendChild(LabelTrack,label,'type',Type{n});
        label  = appendChild(LabelTrack,label,'text',Text{n});
        label  = appendChild(LabelTrack,label,'classentity',classentity{n});
        continuousSamplingPoints=LabelTrack.createElement('continuousSamplingPoints');
        for j =1:3;
            samplePoint = LabelTrack.createElement('samplePoint');
            samplePoint.setAttribute('t',num2str(t(n,j),'%4.3f'));
            samplePoint.setAttribute('y',num2str(y(n,j),'%4.3f'));
            continuousSamplingPoints.appendChild(samplePoint);
        end 
        label.appendChild(continuousSamplingPoints);
        root.appendChild(label);
    end
        xmlwrite(filepath,LabelTrack);   
end
end
function parent  = appendChild(XML,parent,elementname,text)
        element = XML.createElement(elementname);
        element.appendChild(XML.createTextNode(text));
        parent.appendChild(element);
end

% 
%         xmlobj = xmlread(xmlfilename);
%         
%         
%         
%         name            = char(xmlobj.getDocumentElement.getAttribute('name'));
%         externalchange  = char(xmlobj.getDocumentElement.getAttribute('externalchange'));
%         classname       = char(xmlobj.getDocumentElement.getAttribute('classname'));
%         isContinuous    = char(xmlobj.getDocumentElement.getAttribute('isContinuous'));
%         interpolationType = char(xmlobj.getDocumentElement.getAttribute('interpolationType')); 
% 
% 
%         LabelTrack = com.mathworks.xml.XMLUtils.createDocument('LabelTrack');
% 
%             root = LabelTrack.getDocumentElement;
%             root.setAttribute('name',name);
%             root.setAttribute('classname',classname);
%             root.setAttribute('externalchange',num2str(externalchange));
%             
%             root.setAttribute('isContinuous',isContinuous);
%             root.setAttribute('interpolationType',interpolationType);
% 
%             for i=1:size(labels,2)
%                 label = LabelTrack.createElement('label');
% 
%                 startTime=LabelTrack.createElement('starttime');
%                 startTime.appendChild(LabelTrack.createTextNode(num2str(labels(i).startTime)));
%                 label.appendChild(startTime);
%                 endTime=LabelTrack.createElement('endtime');
%                 endTime.appendChild(LabelTrack.createTextNode(num2str(labels(i).endTime)));
%                 label.appendChild(endTime);
%                 timestamp=LabelTrack.createElement('timestamp');
%                 timestamp.appendChild(LabelTrack.createTextNode(num2str(labels(i).timeStamp)));
%                 label.appendChild(timestamp);
%                 value=LabelTrack.createElement('value');
%                 value.appendChild(LabelTrack.createTextNode(num2str(labels(i).value)));
%                 label.appendChild(value);
%                 comment=LabelTrack.createElement('comment');
%                 comment.appendChild(LabelTrack.createTextNode(labels(i).comment));
%                 label.appendChild(comment);
%                 type=LabelTrack.createElement('type');
%                 type.appendChild(LabelTrack.createTextNode(labels(i).type));
%                 label.appendChild(type);
%                 text=LabelTrack.createElement('text');
%                 text.appendChild(LabelTrack.createTextNode(labels(i).text));
%                 label.appendChild(text);
%                 classentity=LabelTrack.createElement('classentity');
%                 classentity.appendChild(LabelTrack.createTextNode(labels(i).classentity));
%                 label.appendChild(classentity);
%                 
%                 continuousSamplingPoints=LabelTrack.createElement('continuousSamplingPoints');
%                 for j =1:3;
%                     samplePoint = LabelTrack.createElement('samplePoint');
%                     samplePoint.setAttribute('t',num2str(labels(i).samplePoint(j,1),'%4.3f'));
%                     samplePoint.setAttribute('y',num2str(labels(i).samplePoint(j,2),'%4.3f'));
%                     continuousSamplingPoints.appendChild(samplePoint);
%                 end 
%                 label.appendChild(continuousSamplingPoints);
% 
%                 root.appendChild(label);
%             end
%             xmlwrite(xmlfilename,LabelTrack);
% 
% output = LabelTrack;