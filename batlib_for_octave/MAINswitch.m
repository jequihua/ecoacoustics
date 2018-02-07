function MAINswitch (data,mode)
switch mode
    case 'all'
        mode = 'ADFSCP';
end
%%
if any (mode == 'A');
    data = getaudioinfo(data);
    savedata(data)
end
%% SONOGRAM CALCULATION
if any (mode == 'D');
    %Energy signal
    data.detdata.WindowTime0      = 2.5e-3;%sec
    data.detdata.Overlap          = 3/4;
    data.detdata.FrequencyFilter  = [12.5 250]*1e3;%Hz
    data.detdata.WindowFunction   = 'hamming';
    data.detdata.MethodKey        = 2;
%         % Energy-based method        
%         data.detdata.Emethod.Threshold      = [12   6];%dB
%         data.detdata.Emethod.PeakWidth      = [50 75];% 0%-100%
%         data.detdata.Emethod.SmoothWindow   = [2.5e-3 2.5e3];% [s Hz]
        % Voice Activity Detection method    
        data.detdata.VADmethod.Threshold    = 0.4;
        % Save mode
        data.detdata.SaveMode = 'd'; %'s': spectrogram,'d' only detections,'s+d' spectrograms and detections;    
    data = detect(data);
    savedata(data);
end
 %%  FEATURE EXTRACTION
 if any (mode == 'F');
    % Feature Extraction
    %Power Spectrum Density
    data.featext.PSD.Window       = 10e-3;%Hz
    data.featext.PSD.Overlap      = 1/2;
    data.featext.PSD.WindowType   = 'rectwin';
    
    %Spectrogram frame
    data.featext.STFT.WindowTime0       = 2e-3;%sec
    data.featext.STFT.Overlap           = 1-1/32;
    data.featext.STFT.WindowType        = 'hamming';
    data.featext.STFT.deltaF            = 70e3;%Hz
    data.featext.STFT.deltaT            = 50e-3;%sec
    
    %Enery-Frequency-Time curve
    data.featext.EFT.Threshold       = 12;%dB
    data.featext.EFT.MaximumArea     = 200;%secHz
    data.featext.EFT.MinimumArea     = 3;%secHz
    data.featext.EFT.AdjustmentStep  = 2;%dB
    data.featext.EFT.FollowFreqPeak  = 0;
    data.featext.EFT.Connectivity    = 4;
    data.featext.EFT.AreaType        = 'Filled';
  
    data.featext.SaveMode ='d';
    data = featext(data);
    savedata(data);
%     %printMain
 end
 %% Crickect Filter;
 if any (mode == 'C');
    load('CricketDet');
    data = CricketFIlter(data,CricketDet);

      % DETECTION SELECTION
    %load(fullfile(pwd,'Basedata','JRS 2015 I 12HRS'));
    data.detdata.Selection.Mode            = 'Area';%% 'MaxE' or 'BestSequence'    
    data.detdata.Selection.AmountLim       = [0 Inf];
    data.detdata.Selection.ClippingFilter  = 0;
    data.detdata.Selection.CricketFilter   = 0;
    data.detdata.Selection.TimeLim         = [0 0];% sec [Ti Tend-Tf]
    data.detdata.Selection.CallDurLim      = [2 40]*1e-3;%sec    
    data.detdata.Selection.AreaLim         = [5 inf]; % Hz*s
    data.detdata.Selection.FiltNaNFeatures = 0; 
    data.detdata.Selection.DetectionLRLim  = [1 inf]; % likelihhod ratio [0 inf]
    data.detdata.Selection.FrequencyFilter = [15 130]*1e3; %Hz;    
    data = selectdetections(data);
    savedata(data)

%% Classify by detection data

%     load('RF_R1_Gen_D')
%     data = classify(data,RF);
%     savedata(data);
        load('RF_R1_Sp_F')
        data = classify(data,RF);
         savedata(data);
end
%% Print worksheets
if any (mode == 'P');
%printMAINaudio
%printMAINcalls
compiletemplate(data)
end
end