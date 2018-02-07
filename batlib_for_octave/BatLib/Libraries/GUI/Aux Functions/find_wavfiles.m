function [wavfiles,directory_name]= find_wavfiles
%This program find all the wav files within a directory includings
%subdirectories and files contained in the subdirectories

[status,mfilepath,messageid] = fileattrib;
path_mfile = mfilepath.Name;

directory_name = uigetdir;
cd (directory_name)
comand = 'dir *.wav/s /b';
[s, pathfiles]= dos(comand);
pathfiles = regexp(pathfiles, '\n','split');
pathfiles(length(pathfiles))=[];

ie=[];
for i = 1:length (pathfiles)
 if isempty(wavread(pathfiles{i},'size'))
     ie = [ie i];
 end
end
pathfiles(ie)=[];

for i = 1:length (pathfiles) 
siz = wavread(pathfiles{i},'size');     
[y, fs, nbits, info] = wavread(pathfiles{i},1);
  
wavfiles(i).name = char(strtok(regexp(pathfiles{i}, '[^\\]*[wavWAV]$', 'match'),'.wav'));
wavfiles(i).path = char(pathfiles{i});
wavfiles(i).samples = siz(1);
wavfiles(i).channels= siz(2);
wavfiles(i).sampling_frequency = fs;
wavfiles(i).bits = nbits;
wavfiles(i).additional_info = info;

end

cd(path_mfile)
end

