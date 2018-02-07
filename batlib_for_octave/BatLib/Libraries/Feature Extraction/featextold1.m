function Data = featextx(Data,inFE,plotkey)
   if nargin == 2
       plotkey=0;
   end

   Data = localfeatext(Data,inFE,plotkey);
   Data = globalfeatext(Data,inFE);
   Data = recountdata (Data);
   savedata(Data);
   
   
   
 