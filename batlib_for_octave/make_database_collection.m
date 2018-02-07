% MAIN INFORMATION
    data.path.Workspace = pwd;    
    data.Title    = ['VAD_test']; 
    data.path.Main = fullfile(data.path.Workspace);
    data.path.Libraries = fullfile(data.path.Workspace,'BatLib');
    data.path.SaveKey  = 1;
    addpath(genpath(data.path.Libraries));

% AUDIO INFORMATION

    data.path.Audio = fullfile(data.path.Workspace,'Audio');
    data.path.WinHeader = data.path.Main;
    data.path.LinHeader = data.path.Workspace;
	data = inidata(data);
	data = getaudioinfo(data);
	savedata(data)
