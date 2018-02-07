function [D,LR,PR,logl] = VADdetector(X,N,vad,ne)
% The algorithm operation is controlled by a small number of parameters:
%
%        pp.of          % overlap factor = (fft length)/(frame increment) [2]
%        pp.ti          % desired output frame increment (10 ms)
%        pp.tj          % internal frame increment (10 ms)
%        pp.ri          % set to 1 to round tj to the nearest power of 2 samples [0]
%        pp.ta          % Time const for smoothing SNR estimate [0.396 seconds]
%        pp.gx          % maximum posterior SNR as a power ratio [1000 = 30dB]
%        pp.gz          % minimum posterior SNR as a power ratio [0.0001 = -40dB]
%        pp.xn          % minimum prior SNR [0]
%        pp.pr          % Speech probability threshold [0.7]
%        pp.ts          % mean talkspurt length (100 ms)
%        pp.tn          % mean silence length (50 ms)
%        pp.ne          % noise estimation: 0=min statistics, 1=MMSE [0]
%
% In addition it is possible to specify parameters for the noise estimation algorithm
% which implements reference [3] from which equation numbers are given in parentheses.
% They are as follows:
%
%        pp.taca      % (11): smoothing time constant for alpha_c [0.0449 seconds]
%        pp.tamax     % (3): max smoothing time constant [0.392 seconds]
%        pp.taminh    % (3): min smoothing time constant (upper limit) [0.0133 seconds]
%        pp.tpfall    % (12): time constant for P to fall [0.064 seconds]
%        pp.tbmax     % (20): max smoothing time constant [0.0717 seconds]
%        pp.qeqmin    % (23): minimum value of Qeq [2]
%        pp.qeqmax    % max value of Qeq per frame [14]
%        pp.av        % (23)+13 lines: fudge factor for bc calculation  [2.12]
%        pp.td        % time to take minimum over [1.536 seconds]
%        pp.nu        % number of subwindows to use [3]
%        pp.qith      % Q-inverse thresholds to select maximum noise slope [0.03 0.05 0.06 Inf ]
%        pp.nsmdb     % corresponding noise slope thresholds in dB/second   [47 31.4 15.7 4.1]
% If convenient, you can call vadsohn in chunks of arbitrary size. Thus the following are equivalent:
%
%                   (a) X=vadsohn(s,fs);
%
%                   (b) [y1,z]=vadsohn(s(1:1000),fs);
%                       [y2,z]=vadsohn(s(1001:2000),z);
%                       y3=vadsohn(s(2001:end),z);
%                       y=[y1; y2; y3];
%
% Note that in all cases the number of output samples will be less than the number of input samples if
% the input length is not an exact number of frames. In addition, if two output arguments
% are specified, the last partial frame samples will be retained for overlap adding
% with the output from the next call to specsub().
%
% Refs:
%    [1] J. Sohn, N. S. Kim, and W. Sung.
%        A statistical model-based voice activity detection.
%        IEEE Signal Processing Lett., 6 (1): 1–3, 1999. doi: 10.1109/97.736233.
%    [2] Ephraim, Y. & Malah, D.
%        Speech enhancement using a minimum-mean square error short-time spectral amplitude estimator
%        IEEE Trans Acoustics Speech and Signal Processing, 32(6):1109-1121, Dec 1984
%    [3] Rainer Martin.
%        Noise power spectral density estimation based on optimal smoothing and minimum statistics.
%        IEEE Trans. Speech and Audio Processing, 9(5):504-512, July 2001.

%      Copyright (C) Mike Brookes 2004
%      Version: $Id: vadsohn.m 3100 2013-06-13 16:05:56Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html


 % default algorithm constants
    dt   = vad.dt;          % true frame increment time
    ta   = vad.ta;
    gmax = vad.gx;          % max posterior SNR = 30 dB
    xn   = vad.xn;          % floor for prior SNR, xi
    gz   = vad.gz;          % floor for posterior SNR
    tn   = vad.tn;
    ts   = vad.ts;
    theta= vad.pr;
    hmm  = vad.hmm;
    
    a=exp(-dt/ta);        % SNR smoothing coefficient
    kk=sqrt(2*pi);        % sqrt(8)*Gamma(1.5) - required constant
    
    
    a01=dt/tn;           % a01=P(silence->speech)
    a00=1-a01;           % a00=P(silence->silence)
    a10=dt/ts;           % a10=P(speech->silence)
    a11=1-a10;           % a11=P(speech->speech)
    
    P01=a10/a01;         % P(silence)/P(speech)
    logP01=log(P01);
%     b11=a11/a10;
%     b01=a01/a00;
%     b00=a01-a00*a11/a10;
%     b10=a11-a10*a01/a00;
    
   
        
    [I,J] = size(X);
    
        % switch vad.ne
        %     case 1
        %     [dp]=estnoiseg(X',dt,ne);	% estimate the noise using MMSE
        %     case 0
        %     [dp]=estnoisem(X',dt,ne);	% estimate the noise using minimum statistics
        % end
        
        switch vad.ge
            case 1                              % Decision Directed
            xu      = 1;                        % dummy unsmoothed SNR from previous frame [2](53)++
            
            case 0                              % Maximum Likelihood
        end
        
        gamma   = max(min(X./N,gmax),gz);  % gamma = posterior SNR [2](10)
        logl    = zeros(I,J);
        logL    = zeros(1,J);
        PR      = zeros(1,J);                % create space for prob ratio vector
        LR      = zeros(1,J);
        for j=1:J                           % loop for each frame
        gammak=gamma(:,j);                  % gamma(j) = a posteriori SNR [2](10)
        
        switch vad.ge
            case 1                          % Decision Directed
                
            xik=a*xu+(1-a)*max(gammak-1,xn);        % xi = smoothed a priori SNR [2](53)
                v=0.5*xik.*gammak./(1+xik);
                gi=(0.277+2*v)./gammak;             % accurate to 0.02 dB for v>0.5
                mv=v<0.5;
                if any(mv)                          % only calculate Bessel functions for v<0.5
                    vmv=v(mv);
                    gi(mv)=kk*sqrt(vmv).*((0.5+vmv).*besseli(0,vmv)+vmv.*besseli(1,vmv))./(gammak(mv).*exp(vmv)); % [2](7)
                end
                xu=gammak.*gi.^2;                   % unsmoothed prior SNR % [2](53)
                
            case 0                          % Maximum Likelihood
            xik = gammak-1;                         
        end
        
        %Likelihood Ratio calculation
        
        logl(:,j) = (xik.*gammak./(1+xik)) - log(1+xik); % log likelihood ratio for each frequency bin [1](3)
        logL(j)  =  sum(logl(:,j))/I;                    % mean log LR -eliminated(over frequency omitting DC term)- [1](4)

        if hmm            
            if j==1
            Gamma = -logP01 +logL(j);            
            else
            Gamma = (a01+a11*Gamma)/(a00+a10*Gamma)*exp(logL(j));            
            end
            logLn = logP01 +log(Gamma);
            LR(j)    = logLn; 
            PR(j) = 2./(1+exp(-logLn  ))-1;            
        else
            LR(j) = logL(j);
            PR(j) = 2./(1+exp(-logL(j)))-1;            
        end        
            %Original code
%             if lggami<0                         % avoid overflow in calculating [1](11)
%                 lggami=lprat+lami+log(b11+b00/(a00+a10*exp(lggami)));
%             else
%                 lggami=lprat+lami+log(b01+b10/(a10+a00*exp(-lggami)));
%             end
%             
%             if lggami<0
%                 
%                 gg=exp(lggami);
%                 prsp(j)=gg/(1+gg);
%             else
%                 prsp(j)=1/(1+exp(-lggami));
%             end
%       
        end
        
    
