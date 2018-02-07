function psdp  = psdparameters(data,psdp,fs,det)

if nargin == 1
    psdp.wt       = data.featext.PSD.Window;
    psdp.ovp      = data.featext.PSD.Overlap;
    psdp.wtype    = data.featext.PSD.WindowType;
    eval(['psdp.whandle = @' psdp.wtype ';'])
elseif nargin == 4
    wn          = nearest2n (psdp.wt,fs);
    psdp.ov      = round (wn*psdp.ovp);
    psdp.w       = window(psdp.whandle,wn);   
    psdp.nfft    = 2*wn;
    psdp.flim    = det.dFp;
end