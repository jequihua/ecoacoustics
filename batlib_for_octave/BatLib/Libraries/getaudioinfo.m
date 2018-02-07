function data = getaudioinfo(data)

    if ispc
        %path0 = pwd;
        cd (data.path.Audio)
        [~, pathfiles]= dos('dir *.wav/s /b /N');
        cd(data.path.Main)
        %cd(path0);
        pathfiles = regexp(pathfiles, '\n','split');
        pathfiles(length(pathfiles))=[];        
    elseif isunix
        folder = changepath(data.path.Audio,data.path);
        if isdir(folder)
            folder = strrep(folder,' ','\ ');
           [~,cmdout] = unix(['find ' folder ' -iname "*.wav"']);
           pathfiles = cellstr(strsplit(cmdout,'\n'));
           if isempty(pathfiles{end})
               pathfiles(end) = [];
           end
        end
    end
    Nf = length(pathfiles);    
    
    Name              =cell(Nf,1);
    Path              =cell(Nf,1);
    Samples           =nan(Nf,1);
    TimeLength        =nan(Nf,1);
    Channels          =nan(Nf,1);
    SamplingFrequency =nan(Nf,1);
    Profundity        =nan(Nf,1);
    AddInfo           =cell(Nf,1);
    Key               =false(Nf,1);
    Error             =cell(Nf,1);
    FileIndex         =zeros(Nf,1);
    

for i = 1:Nf
    if ispc
        strcell = strsplit(pathfiles{i},'\');
    elseif isunix
        strcell = strsplit(pathfiles{i},'/');
    end
    str     = strcell{end};
    iext    = strfind(lower(str),'.wav')-1; 
    Name{i} = char(str(1:iext));
    Path{i} = char(pathfiles{i}(1:end-1));
    disp(['Reading file info ' num2str(i) ' of ' num2str(Nf) ' : ' Name{i}])    
    try    
        afile = Path{i};
        AI = audioinfo(afile);
        
        Samples(i)         = AI.TotalSamples;
        TimeLength(i)      = AI.TotalSamples/AI.SampleRate;
        Channels(i)        = AI.NumChannels;
        SamplingFrequency(i) = AI.SampleRate;
        Profundity(i)      = AI.BitsPerSample;
        AddInfo{i}         = AI;
        Key(i) = 1;
    catch err
        Error{i} = err.message;
        disp(['Error reading file:' Name{i}])
        disp(err.message)
    end
end
    
    data.audioinfo.MainPath          = data.path.Audio;
    data.audioinfo.FilesAmount       = Nf;
    data.audioinfo.Name              = Name;
    data.audioinfo.Path              = changepath(Path,data.path,computer);
    data.audioinfo.Samples           = Samples;
    data.audioinfo.TimeLength        = TimeLength;
    data.audioinfo.Channels          = Channels;
    data.audioinfo.SamplingFrequency = SamplingFrequency;
    data.audioinfo.Profundity        = Profundity;
    data.audioinfo.AddInfo           = AddInfo;
    data.audioinfo.Key               = Key;
    data.audioinfo.Error             = Error;
    data.audioinfo.FileIndex         = FileIndex;
    data.audioinfo.IndexedFilesAmount = nnz(Key);
    data.audioinfo.FileIndex(Key)     = 1:data.audioinfo.IndexedFilesAmount;    
    

        
        

