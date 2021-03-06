function BatLogStruct = BatLogXMLread(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
   xmlobj = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

ElementName = xmlobj.getDocumentElement.getNodeName;
if ~strcmp(ElementName,'BatRecord')
    error('BatLogger signature is not founded in the XML file %s.',filename) 
end
    
   BatLogScheme = loadBLstrcut();
   BatLogStruct = fillBLstruct(xmlobj,BatLogScheme);

function BatLogScheme = loadBLstrcut()
GPS         = struct(   'Valid',        @(x) char(x),...
                        'Position',     @(x) sscanf(char(x),'%f %f'),...
                        'Altitude',     @(x) sscanf(char(x),'%f'),...
                        'CH1903',       @(x) char(x),...
                        'HDOP',         @(x) sscanf(char(x),'%f'),...
                        'SatsUsed',     @(x) int8(sscanf(char(x),'%d')),...
                        'GPSTimestamp', @(x) char(x),...
                        'GPSAge',       @(x) char(x)...
                        );
                    
Trigger     = struct(   'TRIG_MODE',        @(x) char(x),...
                        'Version',          @(x) char(x),...
                        'Event',            @(x) char(x),...
                        'PRETRIG_TIME_MS',  @(x) int16(sscanf(char(x),'%d')),...
                        'POSTTRIG_TIME_MS', @(x) int16(sscanf(char(x),'%d')),...
                        'TRIG_PAR0',        @(x) int16(sscanf(char(x),'%d')),...
                        'TRIG_PAR1',        @(x) int16(sscanf(char(x),'%d')),...
                        'TRIG_PAR2',        @(x) int16(sscanf(char(x),'%d')),...
                        'TRIG_PAR3',        @(x) int16(sscanf(char(x),'%d')),...
                        'TRIG_PAR4',        @(x) int16(sscanf(char(x),'%d')),...
                        'TRIG_PAR5',        @(x) int16(sscanf(char(x),'%d')),...
                        'TrigValue0',       @(x) int16(sscanf(char(x),'%d')),...
                        'TrigValue1',       @(x) int16(sscanf(char(x),'%d')),...
                        'TrigValue2',       @(x) int16(sscanf(char(x),'%d')),...
                        'TrigValue3',       @(x) int16(sscanf(char(x),'%d')),...
                        'TrigValue4',       @(x) int16(sscanf(char(x),'%d')),...
                        'TrigValue5',       @(x) int16(sscanf(char(x),'%d'))...
                        );
                    
BatLogScheme = struct(  'Firmware',    @(x) char(x),...
                        'SN',          @(x) int16(sscanf(char(x),'%u')),...
                        'Filename',    @(x) char(x),...
                        'DateTime',    @(x) datenum(char(x),'dd.mm.yyyy HH:MM:SS'),...
                        'Duration',    @(x) sscanf(char(x),'%f'),...
                        'Samplerate',  @(x) sscanf(char(x),'%f'),...
                        'Temperature', @(x) int8(sscanf(char(x),'%f')),...
                        'BattVoltage', @(x) single(sscanf(char(x),'%f')),...
                        'GPS',         GPS,...
                        'Trigger',     Trigger...
                        );
function scheme = fillBLstruct(xmlobj,scheme)
fields = fieldnames(scheme);
I = length(fields);
    for i = 1:I
        if ~isstruct(scheme.(fields{i}))
        element = xmlobj.getElementsByTagName(fields{i});
           if element.getLength
                obj    = element.item(0).item(0).getData;
                fcn    = scheme.(fields{i});
                data   = fcn(obj);
                scheme.(fields{i}) = data;
           else
                scheme.(fields{i}) = '';
                %disp(['Unable to convert data from Node: ' fields{i}])
           end        
        else
            subscheme = fillBLstruct(xmlobj,scheme.(fields{i}));
            scheme.(fields{i}) = subscheme;
        end
    end

