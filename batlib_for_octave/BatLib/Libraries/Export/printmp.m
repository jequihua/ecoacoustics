classdef printmp
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pageSize
        pageH
        pageL
        
        frameH
        frameL
        framePos
        
        track
        trackN
        trackH
        trackL
        trackPos
        
        trackCmap
        trackClim
        
        trackZlim
        
        trackTsize
        trackTlim
        trackTtick
        trackTtickstr
            
        trackFsize   
        trackFlim
        trackFtick
        trackFtickstr
        
        pageTitle
        pageNote1
        pageNote2
        pageNumber
        pageNfig
        
        savename
           
            
        
%         frameL
%         P.Psiz     = inPrint.Page.Size;
%         P.H        = inPrint.Page.Size(1);
%         P.L        = inPrint.Page.Size(2);
%         P.h        = inPrint.Page.FrameHeigth;
%         P.l        = inPrint.Page.FrameLength;
%         P.pos      = [ P.l, P.h, P.L-2*P.l, P.H-2*P.h];
%         P.ln.N   = inPrint.PlotA.LinesxPage;
%         P.ln.H   = P.pos(4)/P.Nline*0.85;
%         P.ln.L   = P.pos(3);
%         P.ln.pos = @(l) [P.l,(P.H-P.h-l*P.Hline), P.Lline, P.Hline];  
%         P.ln.dt  = inPrint.PlotA.SpecTimeLength;
%         %Pg.ln.t(2) = inPrint.PlotA.DetTimeLength;
%         P.ln.tlim  = @(l) [(l-1) l]*P.ln.dt;
%         P.ln.ttick = 0:0.05:P.ln.dt;
%         P.ln.ttickstr = @(tlim) tlim(1):0.05:tlim(2);
% 
%         P.ln.flim  = @(fs)[0 min(fs/2,150e3)];
%         P.ln.ftick = @(flim)flim(1):20e3:flim(2);
%         P.ln.ftick = P.ln.flim/1000;
% 
% 
%         P.text.title  =  @(n,name) strcat('File ',num2str(n),': ',name);
%         P.text.note1  =  @(fs,wt,ov,df,dt,N) strcat('Fs = ',num2str(fs/1e3,'%2.1f kHz'),...
%                                                   ', wt = ',num2str(wt*1e3,'%2.1f ms'),...
%                                                   ', ov = ',num2str(ov/wt*100,'%%%2.2f'),...
%                                                   ', df = ',num2str(df,'%2.1 Hz'),...                                              
%                                                   ', dt = ',num2str(dt*1e3,'%2.1 ms'),...
%                                                   'Ndet = ',num2str(N));
%         P.text.page    = @(p) num2str(p);
%         P.text.nfile   = @(n,N) strcat(num2str(n),'/',num2str(N));
%         P.text.pdfname = inPrint.Filename;
    end
    
    methods 
        function obj = printmp(pgsiz,frsiz,N)
           
            H = pgsiz(1);
            L = pgsiz(2);
            h = frsiz(1);
            l = frsiz(2);            
            obj.pageSize = pgsiz;
            obj.pageH    = H;
            obj.pageL    = L;            
            obj.frameH   = h;            
            obj.frameL   = l;
            obj.framePos = [l h L-2*l H-2*h];            
            obj.trackN   = N;            
        end
        
        function settrack(obj,N)
        end
        
        
    end
    
end

