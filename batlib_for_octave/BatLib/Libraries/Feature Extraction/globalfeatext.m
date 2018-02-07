function Data = globalfeatext(Data,inFE)
 % Features Name    
    Data.FeatListStr = [Data.FeatListStr Gfeatures()'];
    Data.FeatList    = 1:size(Data.FeatListStr,2);
 % Feature Extraction
    X = Gfeatures ([],Data);    
    Data.Data = [Data.Data X];
    
 
