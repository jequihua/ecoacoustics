function scalartrack = data2strackRAW(data,nfile,scalartrack)

    filepath =  data.Atlas.ProjectPath{nfile}.datatracks;
    Ntrack = length(scalartrack);
 
    for n = 1:Ntrack
        switch scalartrack{n}.Attributes.name
        
            case 'DetSignal'
                
                D          = data.detdata.VADmethod.DetectionSignal{nfile};
                J          = size(D,2);
                dT         = data.detdata.TimeStep(nfile)*1e3;
                Tmin       = data.detdata.TimeZero(nfile)*1e3; 
                T          = (0:J-1)*dT + Tmin;
                
                X            = nan(2,J); 
                X(1,1:J)     = T;
                X(2,1:J)     = D; 
                
                scalartrack{n}.Attributes.min          = num2str(min(D),'%10.8f');
                scalartrack{n}.Attributes.max          = num2str(max(D),'%10.8f');
                
            case ''
                
                
        end        
        
        filename = fullfile(filepath,scalartrack{n}.Text);
        writeraw(X,filename)
        %save(fullfile(filepath,'DetSignal.mat'),'X')                
    end
end
 
function writeraw(X,filename)
precision   = 'double'; 
skipdata    = 0;
machinefmt  = 'b';

fileID = fopen(filename,'w+');
fwrite(fileID,X,precision,skipdata,machinefmt);
fclose(fileID);
end