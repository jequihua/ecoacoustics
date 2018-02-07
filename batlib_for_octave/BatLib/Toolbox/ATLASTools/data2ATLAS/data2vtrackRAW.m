function vectrack = data2vtrackRAW(data,nfile,vectrack)

    filepath =  data.Atlas.ProjectPath{nfile}.datatracks;
    Ntrack = size(vectrack,2);
 
    for n = 1:Ntrack
        switch vectrack{n}.Attributes.name
        
            case 'Spectrogram'
                SpecData   = load(data.detdata.SavePath{nfile},'p','t','f');
                P0         = SpecData.p;
                T0         = SpecData.t*1e3;
                F0         = SpecData.f/1e3;
                % Reshape                
                  rho           = 1000/500;% pixels per miliseconds (pix/ms);
                  npixels       = round((T0(end)-T0(1))*rho);
                  T             = linspace(T0(1),T0(end),npixels);
                  
                  mpixels       = str2double(vectrack{n}.Attributes.trackHeight);
                  Fmin          = F0(1);
                  Fmax          = F0(end);
                  F             = linspace(Fmin,Fmax,mpixels);
                  P             = interp2(T0',F0,P0,T,F','linear');
                  
                PdB          = pow2db(P); 
                [I,J]        = size(PdB);
                X            = nan(I+1,J); 
                X(1,1:J)     = T;
                X(2:I+1,1:J) = PdB(I:-1:1,1:J); 
                
                vectrack{n}.Attributes.min          = num2str(nanmean(PdB(:)),'%10.8f');
                vectrack{n}.Attributes.max          = num2str(max(PdB(:)),'%10.8f');
                vectrack{n}.Attributes.dimension    = num2str(I,'%10.0f');
                
                
            case 'Call Spectrogram'
                
                
        end        
        
        filename = fullfile(filepath,vectrack{n}.Text);
        writeraw(X,filename)
        %save(fullfile(filepath,'Spectrogram.mat'),'X')       
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