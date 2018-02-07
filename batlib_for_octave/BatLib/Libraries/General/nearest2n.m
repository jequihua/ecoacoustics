function w = nearest2n (wt,fs)
        w0 = wt*fs;
        a = floor(log2(wt*fs));
        [~,b] = min([w0-2.^a 2.^(a+1)-w0]);
    w = 2.^(a+b-1);