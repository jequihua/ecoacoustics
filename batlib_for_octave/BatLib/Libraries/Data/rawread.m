function [X,nbits,samples] = rawread(filepath,N, Encoding, Skip, ByteOrder)

    FileExt = regexp(filepath,'.raw$','once');
if isempty(FileExt)
    error('Unknown file extension')
end

switch nargin 
    case 4
        ByteOrder = 'l';
    case 3 
        ByteOrder = 'l';    
        Skip = 0;
    case 2
        ByteOrder = 'l';    
        Skip = 0;
        Encoding = 'int16';
    case 1
        ByteOrder = 'l';    
        Skip = 0;
        Encoding = 'int16';
        N = inf;       
end
    try
    outputformat = [Encoding '=>double'];    
    fileID  = fopen(filepath); 
    X0      = fread(fileID, N, outputformat, Skip, ByteOrder);
    fclose(fileID);
    X0max   = double(intmax(Encoding)); 
    X = X0/X0max;
    
    nbits    =  16;
    samples =  length(X);
    catch
        fclose('all');
    end
end

