function XMLProject = data2ProjectXML( data,nfile,ProjectPath,pdata)
    templatepath = fullfile(data.path.Libraries,'Libraries','ATLASTools','Template','T.xml');
    XMLProject = xmlread(templatepath);
    %XMLProject = com.mathworks.xml.XMLUtils.createDocument('AnnotationProject');
        
    ScreenSize = get(0,'screensize');
    TotalTime  = ceil(data.audioinfo.TimeLength(nfile)*1e3);
    Zoom       = 1;%TotalTime/ScreenSize(3);
    
    name       = data.audioinfo.Name{nfile};
    length     = num2str(TotalTime,'%10.1f');
    zoom       = num2str(Zoom,'%10.8f');
    command    = '';
    
    filepath = fullfile(ProjectPath.main,[name '.xml']);
        
    
	root = XMLProject.getDocumentElement;
    root.setAttribute('name',name);
    root.setAttribute('length',length);
    root.setAttribute('zoom',zoom);
    root.setAttribute('commandlineExecutionString',command);

    vectrack = pdata.vectrack;
    scatrack = pdata.scatrack;
    XMLProject = settrack(XMLProject,vectrack);
    XMLProject = settrack(XMLProject,scatrack);
    xmlwrite(filepath,XMLProject);

end
    
function XML = settrack(XML,track)
    trackname = inputname(2);        
    child  = XML.getElementsByTagName(trackname);
    nchild = child.getLength;
    for i = 0:(nchild-1)
        for j = 1:size(child,2)
            atrname = track{j}.Attributes.name;
           if strcmp(atrname, char(child.item(i).getAttribute('name')));
                field = fieldnames(track{j}.Attributes);
                for k = 1: size(field,1)
                    atrname     = field{k};
                    atrvalue    = track{j}.Attributes.(field{k});
                    XML.getElementsByTagName(trackname).item(i).setAttribute(atrname,atrvalue);
                end
           end
        end
    end
end



    
function parent  = appendChild(XML,parent,elementname,struct)
for i = 1:size(struct,2)    
    if ~isempty(struct{i}.Text)
        
        element = XML.createElement(elementname);
        element.appendChild(XML.createTextNode(struct{i}.Text));
        
        name = fieldnames(struct{i}.Attributes);
        for j = 1:length(name) 
        element.setAttribute(name{j},struct{i}.Attributes.(name{j}));      
        end
        
        parent.appendChild(element);
    end
end
end

