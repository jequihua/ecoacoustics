function [out] = getvadparam (dT,theta)
      VAD.dt  = dT;             % true frame increment time
      VAD.pr  = theta;            % Speech probability threshold
      VAD.gx  = 10^(30/10);     % maximum posterior SNR = 30dB
      VAD.gz  = 10^(-40/10);    % minimum posterior SNR = -40dB

      VAD.ne  = 0;              % noise estimation: 0=min statistics, 1=MMSE [0]

      VAD.ge  = 1;                % xi estimation: 0= Itakura Saito (ISD), 1 = Decision Directed (DD) [1]
      VAD.ta  = -dT/log(0.98);    % Time const for smoothing SNR estimate = -tinc/log(0.98) from [2]
      VAD.xn  = 0;                % minimum prior SNR = -Inf dB

      VAD.hmm = 1;                % HMM-Based Hang Over [0]
      VAD.ts  = 0.04;             % mean talkspurt length (100 ms)
      VAD.tn  = 0.005;             % mean silence length (50 ms)


      % Estimate noise spectrum using minimum statistics
      NE.taca   =-dT/log(0.95);       % smoothing time constant for alpha_c = -tinc/log(0.7) in equ (11)
      NE.tamax  =-dT/log(0.96);      % max smoothing time constant in (3) = -tinc/log(0.96)
      NE.taminh =-dT/log(0.3);       % min smoothing time constant (upper limit) in (3) = -tinc/log(0.3)
      NE.tpfall =0.064/10;           % time constant for P to fall (12)
      NE.tbmax  =-dT/log(0.8);             % max smoothing time constant in (20) = -tinc/log(0.8)
      NE.qeqmin =2;                  % minimum value of Qeq (23)
      NE.qeqmax =14;                 % max value of Qeq per frame
      NE.av     =2.12;               % fudge factor for bc calculation (23 + 13 liNEs)
      NE.td     =0.01;              % time to take minimum over
      NE.nu     =16;                  % number of subwindows
      NE.qith   =[0.03 0.05 0.06 Inf];% noise slope thresholds in dB/s
      NE.nsmdb  =[47 31.4 15.7 4.1];

      out.vad = VAD;
      out.ne  = NE;