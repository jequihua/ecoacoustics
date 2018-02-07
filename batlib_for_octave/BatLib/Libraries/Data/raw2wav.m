function raw2wav(path)

info = dir(path);
I = size(info,1);
fs = 500e3;

for i = 1:I
    if ~info(i).isdir;
    filename = info(i).name;
    n = length(filename);    
        if strcmp(filename(n-3:n),'.raw')
            disp(['Converting file: ' filename]);
            try
            filepath = fullfile(path,filename);
            [X,nbits,~] = rawread(filepath);
            %plot(X);
            wavfilepath = strrep(filepath,'.raw','.wav');
            
                audiowrite(wavfilepath,X,fs,'BitsPerSample',nbits);
                
            catch err
                disp(err.message)
            end
        end
    end
end




