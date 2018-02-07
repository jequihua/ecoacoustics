function det = detparameters(data,i)

    det.flim       = data.detdata.FrequencyFilter;
    det.wt         = data.detdata.WindowTime(i);
    det.wn         = data.detdata.WindowSize(i);
    det.Smax       = data.audioinfo.Samples(i);
    
    fs             = data.audioinfo.SamplingFrequency(i); 
    det.fs         = fs;
    
    
    switch data.detdata.MethodKey
        case 1
            det.S       = round(data.detdata.Emethod.DetectionTimeInt{i}*fs);
            det.dS      = round(data.detdata.Emethod.DetectionTimeLength{i}*fs);
            det.Tp      = data.detdata.Emethod.DetectionTimePoint{i};
            det.dTp     = data.detdata.Emethod.DetectionTimeInt{i};
            det.dFp     = data.detdata.Emethod.DetectionFreqInt{i};
            
        case 2
            det.S       = round(data.detdata.VADmethod.DetectionTimeInt{i}*fs);
            det.dS      = round(data.detdata.VADmethod.DetectionTimeLength{i}*fs);
            det.Tp      = data.detdata.VADmethod.DetectionTimePoint{i};
            det.dTp     = data.detdata.VADmethod.DetectionTimeInt{i};
            det.dFp     = data.detdata.VADmethod.DetectionFreqInt{i};
            
    end
    
    det.Sp      = round(det.Tp*fs);
    det.dSp     = det.wn;    
    
    
