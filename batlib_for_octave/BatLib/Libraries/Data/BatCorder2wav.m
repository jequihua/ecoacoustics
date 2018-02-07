function BatCorder2wav

    folder = uigetdir('','Select folder with RAW data');
    
    disp('Searching Raw data in folders')
    [~, pathfiles]= dos(['dir ' folder  '\*.raw /s /b /N']);
    pathfiles = regexp(pathfiles, '\n','split');
    pathfiles(length(pathfiles))=[];
    Nf = length(pathfiles);
        
    
    Name              =cell(Nf,1);
    Path              =cell(Nf,1);
for i = 1:Nf
    strcell = strsplit(pathfiles{i},'\');
    str     = strcell{end};
    iext    = strfind(lower(str),'.raw')-1; 
    Name{i} = char(str(1:iext));
    Path{i} = char(pathfiles{i});
end



fs = 500e3;

for i = 1:Nf    
            filename =  Name{i};    
            disp(['Converting file: ' num2str(i) ' of '  num2str(Nf) '-' filename]);
            try
            filepath = Path{i};
            [X,nbits,~] = rawread(filepath);
            %plot(X);
            wavfilepath = strrep(filepath,'.raw','.wav');
            
                audiowrite(wavfilepath,X,fs,'BitsPerSample',nbits);
                
            catch err
                disp(err.message)
            end

end