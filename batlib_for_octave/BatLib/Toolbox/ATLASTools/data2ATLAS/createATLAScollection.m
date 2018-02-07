function XMLCollection = createATLAScollection( data,mode)
if nargin == 1
    mode = 'default';
end
switch mode
    case 'default'
% Routine parameters
    N        = data.audioinfo.FilesAmount;
    if isfield(data,'featext')
        Key  = data.featext.Key;
    elseif isfield(data,'detdata')
        Key  = data.detdata.Key;
    else
        Key  = data.audioinfo.Key;
    end
end
	templatepath = fullfile(data.path.Libraries,'Libraries','ATLASTools','Template','C.xml');
    XMLCollection = xmlread(templatepath);
    root = XMLCollection.getDocumentElement;
    
for n = 1:N;
    if Key(n)
        disp(['Creating AtlasProject: ' num2str(n)]);
        try
        [data,XMLProject] = createATLASproject(data,n);
        name = [char(XMLProject.getDocumentElement.getAttribute('name')) '.xml'];
        path = data.Atlas.ProjectPath{n}.main;
        filename = fullfile(path,name);
        Project = XMLCollection.createElement('Project');
        Project.setAttribute('name',name);
        Project.setTextContent(filename);
        root.appendChild(Project);
        end
    end
end
    filepath = fullfile(data.path.Atlas,'ATLAScollection.xml');
    xmlwrite(filepath,XMLCollection);
end

