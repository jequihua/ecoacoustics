function XMLclass = data2ClassXML(data,lclass,ProjectPath)
    Nctrack         = size(lclass,2);
for i = 1:Nctrack
    filepath = fullfile(ProjectPath.classes,lclass{i}.Text);
    templatepath = fullfile(data.Atlas.TemplatePath,'classes',lclass{i}.Text);
    XMLclass = xmlread(templatepath);
    
    switch lclass{i}.Attributes.name        
        case 'generic'            
%             entity{1}.name.Text     = 'labeled';
%             entity{1}.color.Text    = '-4144960';
%             entity{1}.continuousColor.Text    = '-12566464';
%             entity{1}.id.Text       = '0';
        case 'Species'
            if isfield(data,'classinfo')
                SpeciesList = data.classinfo.SpeciesList;
                else
                SpeciesList = {'undefined';'A';'B';'C'};
            end
                L = size(SpeciesList,1);
                cmap = round(lines(L)*255);
                SpeciesColor = cell(1,L);
                for l = 1:L
                    SpeciesColor{l} = num2str(rgb2cnum(cmap(l,:,:)));
                end
                    continuousColor(1:L) = {'-16777216'};
                    id                   = cellstr(num2str((1:L)'));
                    
                    root = XMLclass.getDocumentElement;
                for l = 1:L
                entity  = XMLclass.createElement('entity');
                entity  = appendChild(XMLclass,entity,'name',SpeciesList{l});
                entity  = appendChild(XMLclass,entity,'color',SpeciesColor{l});
                entity  = appendChild(XMLclass,entity,'continuousColor',continuousColor{l});
                entity  = appendChild(XMLclass,entity,'id',id{l});
                root.appendChild(entity); 
                end
            
            
        case 'Genus'
            
        case 'Familiy'
            
        case 'Call Phase'
            
        case 'Detection'
    end
        xmlwrite(filepath,XMLclass);
end
end

function parent  = appendChild(XML,parent,elementname,text)
        element = XML.createElement(elementname);
        element.appendChild(XML.createTextNode(text));
        parent.appendChild(element);
end
    