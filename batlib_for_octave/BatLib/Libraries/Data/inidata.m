function data = inidata(indata)

    data = indata;
    try
        pathMAIN = changepath(data.path.Main,data.path);
        data.path.Main = pathMAIN;
        disp(['New Main path : ' pathMAIN]); 
        data.path.Audio = changepath(data.path.Audio,data.path);
        disp(['New Audio path : ' data.path.Audio]); 
    catch
        pathMAIN = data.path.Main;
    end
    %data.path.Libraries  = fullfile(data.path.Main,'Libraries');
    data.path.Basedata   = fullfile(pathMAIN,'Basedata');
    data.path.Printouts  = fullfile(pathMAIN,'Printouts');
    data.path.Worksheets = fullfile(pathMAIN,'Worksheets');
    data.path.Sonograms  = fullfile(pathMAIN,'Sonograms');
    data.path.Sytem      = computer;
    %plib = genpath(data.path.Libraries);
    %addpath (plib);
    
    pathname = fieldnames(data.path);
    for  i = 1:length(pathname)
        switch (pathname{i})
            case{'Basedata','Printouts','Worksheets','Sonograms'}
            if ~isdir(data.path.(pathname{i}))
                mkdir(data.path.(pathname{i}))
            end
        end
    end