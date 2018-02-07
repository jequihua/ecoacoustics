function BCstruct = BatCorTXTread(txtpath)
try
   fileID = fopen(txtpath);
   C = textscan( fileID, '%s','delimiter','\n','MultipleDelimsAsOne',1,'CollectOutput',1);
   txt = C{1};
catch
   error('Failed to read TXT file %s.',txtpath);
end

Signature = regexp(txt{1},'^Batcorder');
if isempty(Signature) || Signature~=1
    error('BatCorder signature is not founded in the XML file %s.',txtpath) 
end

BCstruct = fillBCstruct(txt,txtpath);


    function BCstruct = fillBCstruct(txt,txtpath)
    BCstruct = loadBCstruct();

    Header = strsplit(txt{1},' / ');
    BCstruct.Firmware           = sscanf(Header{1},'Batcorder %s');
    BCstruct.LogFile.Version    = sscanf(Header{2},'logfile %s');
    BCstruct.LogFile.Date       = txt{2}(12:end);
    BCstruct.LogFile.FilePath   = txtpath;
    BCstruct.Samplerate         = 500e3;

    iFile = findrow (txt,'^[MAT]\t');
if ~isempty(iFile)
    for j = 1:length(iFile)
        %info = strsplit(txt{iFile(j)},'\t');
        info0 = textscan(txt{iFile(j)},'%s','delimiter','\t');
        info = info0{1};
        BCstruct.RecordingMode{j,1}   = info{1};
        BCstruct.DateTime(j,1)        = datenum([info{2} info{3}],'dd.mm.yyHH:MM:SS');
        BCstruct.Filename{j,1}        = info{4};
        BCstruct.Duration{j,1}        = info{5};
    end
end

    iA  = findrow (txt,'^Auto ');
    iT =  findrow (txt,'^Timer ');

    iMode = sort([iA; iT]);
    k     = 1;
if ~isempty(iMode)
    for j = 1:2:length(iMode)
    infoOn  = strsplit(txt{iMode(j)},'\t');   
    infoOff = strsplit(txt{iMode(j+1)},'\t');
    param = int16(sscanf(infoOn{5},'"%d;%d;%d;%d"'));
    
    BCstruct.Setup.Mode{k,1}              = sscanf(infoOn{1},'%s on');
    BCstruct.Setup.OnOff{k,1}             = [infoOn{2} ' ' infoOn{3} ' - ' infoOff{2} ' ' infoOff{3}];
    BCstruct.Setup.Filecode{k,1}          = infoOn{4};
    BCstruct.Setup.Quality(k,1)           = param(1);
    BCstruct.Setup.Threshold(k,1)         = param(2);
    BCstruct.Setup.Posttrigger(k,1)       = param(3);
    BCstruct.Setup.CriticalFrequency(k,1) = param(4);    
    k=k+1;
    end
end
    


    function BCstruct = loadBCstruct()

Logfile      =  struct( 'Version',      '',...
                        'FilePath',     '',...
                        'Date',         '',...
                        'Notes',        ''...
                        );
Setup       =   struct( 'Mode',         '',...
                        'OnOff',        '',...
                        'Filecode',     '',...
                        'Quality',      [],...
                        'Threshold',    [],...
                        'Posttrigger',  [],...
                        'CriticalFrequency',[]...
                        );
                    
BCstruct    =  struct(  'Firmware',     '',...
                        'LogFile',      Logfile,...
                        'Setup',        Setup,...
                        'Filename',     '',...
                        'RecordingMode','',...
                        'DateTime',     [],...
                        'Duration',     [],...
                        'Samplerate',   []...
                        );
    function index = findrow (txt,exp)
                C = regexp(txt,exp);
                index = find(cellfun(@(x) ~isempty(x),C,'UniformOutput',1));
  