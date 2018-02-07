function [P,F,T] = spec(x,wn,ov,Nfft,fs,wf,mode)
    if nargin ==6
        mode = 'tftb';
    end
    
    switch mode
        case 'matlab'
            if isempty(wf);
                wf = 'hamming';
            end
                w = window(wf,wn);            
        [~,F,T,P] =  spectrogram (x,w,ov,Nfft,fs);
        
        case 'tftb'
            
        if isrow(x)
            x = x';
        end           
            
        Nx = length(x);
        
        if isempty(wf);
                wf = 'hamming';
        end
                w = window(wf,wn+1);                
%         if length(wn)>1;
%             wn = length(wn)- rem(length(wn),2);
%             win = wn;
%         else
%             wn = wn;
%             win =  tftb_window(odd(wn),'Hamming');            
%         end
        
        if isempty(Nfft)
           N = wn;
        else 
           N = 2*Nfft;
        end
        
        n = wn/2+1:wn-ov:Nx-wn/2+1;
        %spectrogram
        [s,n,f] = tfrsp(x,n,N,w,0);
        
        T = (n-1)/fs;
        
        iF  = f>=0;
        F = f(iF)*fs;
        
        a = 2/fs;
        P = a*s(iF,:);
    end
end

