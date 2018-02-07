function [p,f,t,i0,j0,P,F,T] = stftframe(audio,stft,m)
%STFT calculation
        [P,F,T] = spec (audio,stft.wn,stft.ov,[],stft.fs,stft.wf,'tftb');
% Spectrogram Reference point 
        iFb = interval(F,stft.fplim(m,:),'[]');
        iTb = interval(T,stft.tplim(m,:),'[]');
        P0 = P(iFb,iTb);
        i0 = P == max(P0(:)) & logical(single((iFb)*single(iTb)));
        [iF0,iT0] = find(i0,1,'first');
        T0 = T(iT0);
        F0 = F(iF0);        
%         [maxE it0] = max(sum(P(iFb,iTb),1),[],2);
%         iT0 = it0 + iTb(1)-1;      
%         
%         [P0 if0] = max(P(iFb,iT0,[],1));
%         iF0 = if0 + iFb(1)-1;
      
% EFT surface delimitation
        iF = interval(F,F0 +[-1 1]*stft.df/2,'[]');
        iT = interval(T,T0 +[-1 1]*stft.dt/2,'[]');
        p = P(iF,iT);
        f = F(iF);
        t = T(iT);        
        
        i0 = find(f==F0);
        j0 = find(t==T0);
        
       