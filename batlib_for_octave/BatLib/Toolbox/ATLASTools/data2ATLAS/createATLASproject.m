function [data,XMLProject] = createATLASproject(data,nfile)

% Create folders
    if ~isfield(data.path,'Atlas')
        data.path.Atlas = fullfile(data.path.Main,'Atlas');
        plib = genpath(data.path.Libraries);
        addpath (plib);
        if ~isdir(data.path.Atlas)
            mkdir(data.path.Atlas)
        end
        N = data.audioinfo.FilesAmount;
        data.Atlas.ProjectPath = cell(N,1);
    end
        ProjectName             = data.audioinfo.Name{nfile};%num2str(nfile);%
        ProjectPath.main        = fullfile(data.path.Atlas,ProjectName);
        ProjectPath.classes     = fullfile(ProjectPath.main,'classes');
        ProjectPath.datatracks  = fullfile(ProjectPath.main,'datatracks');
        ProjectPath.labeltracks = fullfile(ProjectPath.main,'labeltracks');
        ProjectPath.media       = fullfile(ProjectPath.main,'media');
         
        folders = fieldnames(ProjectPath);
        for i = 1:size(folders,1)
            if ~isdir(ProjectPath.(folders{i}))
                mkdir(ProjectPath.(folders{i}));
            end
        end
            
        data.Atlas.ProjectPath{nfile}  = ProjectPath;
        data.Atlas.TemplatePath        = fullfile(data.path.Libraries,'Libraries','ATLASTools','Template');
        
% Load data and replace
    % vector track
    vectrack{1}.Attributes.name = 'Spectrogram';
    vectrack{1}.Text            = 'Spectrogram.raw';
    vectrack{1}.Attributes.trackHeight = '200';
    pdata.vectrack = data2vtrackRAW(data,nfile,vectrack);
    % scalar track
    scatrack{1}.Attributes.name = 'DetSignal';
    scatrack{1}.Text            = 'DetSignal.raw';
    pdata.scatrack = data2strackRAW(data,nfile,scatrack);

% Create classes
    lclass{1}.Attributes.name               = 'generic';
    lclass{1}.Text                          = 'generic.xml';
    lclass{2}.Attributes.name               = 'Species';
    lclass{2}.Text                          = 'Species.xml';
    XMLclass = data2ClassXML(data,lclass,ProjectPath);   
    
% Create  label track
    ltrack{1}.Text               = 'SpID.xml';
    XMLLabel = data2LabelXML (data,nfile,ltrack,ProjectPath);
    
% Modify XML project
    XMLProject = data2ProjectXML( data,nfile,ProjectPath,pdata);
    


    
