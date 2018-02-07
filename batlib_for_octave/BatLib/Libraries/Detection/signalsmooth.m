function y = signalsmooth(x,window,deltaT)
            w  = round(window/deltaT); 
            y = smooth(x,w);