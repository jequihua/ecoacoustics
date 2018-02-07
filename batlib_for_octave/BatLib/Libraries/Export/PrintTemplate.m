classdef PrintTemplate < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        pagehandle        
        pageSize = [21 29.7];
        pageH    = 21;
        pageL    = 29.7;
        pageTitle = '';
        pageNote1 = '';
        pageNote2 = '';
        pageNumber = 0;
        
        
        frameH   = 21-1.5;
        frameL   = 29.7-1.5;
        framel   
        frameh   
        framePos        
        
        trackhandle
        trackNumber
        trackNtotal
        trackH
        trackL
        trackPos
        trackTitleColor
        
            trackCmap
            trackClim
            trackClimMode
            trackZlim

            trackTsize
            trackTlim        
            trackTtick
            trackTtickdist
            trackTtickstr

            trackFsize   
            trackFlim
            trackFtick
            trackFtickdist
            trackFtickstr
            
        boxcolor
            
        audiofileName
        audiofileNumber
        audiofileNtotal
        
        savename
        savekey = false;

    end
    
    methods 
        function tmp = PrintTemplate(page,frame,N)
           
            if nargin ==0                
                N    = 4; 
                page = [29.7 21];
                frame= page - 2*[1.5 1.5];
            elseif nargin==1
                N =page;
                page = [29.7 21];
                frame= page - 2*[1.5 1.5];
            end
                L = page(1);
                H = page(2);            
                l = frame(1);
                h = frame(2);
            
            tmp.pageSize = [L,H];
            tmp.pageH    = H;
            tmp.pageL    = L;            
            tmp.frameH   = h;            
            tmp.frameL   = l;
            tmp.frameh   = (H-h)/2;
            tmp.framel   = (L-l)/2;
            tmp.framePos = [(L-l)/2 (H-h)/2 l h];            
            
            tmp.trackNumber   = 1;
            tmp.trackNtotal   = N;
            tmp.trackH        = h/N;
            tmp.trackL        = l;
            tmp.trackPos      = zeros(N,4);
            tmp.trackCmap     = 1-gray(128);
            tmp.trackClimMode = 'mean';
            for n = 1:N
               tmp.trackPos(n,:) = [tmp.framePos(1), tmp.framePos(2)+(N-n)*tmp.trackH,...
                                    tmp.trackL     , tmp.trackH*0.85];
            end
            
        end
        function  h = newfig(tmp,h)
            if nargin <2
                h = figure;
            else
                figure(h);
            end
            
            set(h,'Units', 'centimeters','Position',[2 2 tmp.pageSize])
            tmp.pagehandle = h;        
            tmp.newtracks(h);
        end
        
        function newtracks(tmp,h)
            if nargin <2
                h = tmp.pagehandle;
            end
            figure(h);
            N = tmp.trackNtotal;
            for n = 1:N
               %t   = subplot(N,1,n);
               t = axes('Units','centimeters','Position',tmp.trackPos(n,:)); 
               tmp.trackhandle(n)=t;
            end            
        end
        
        function set.trackTsize(tmp,dt)
            tmp.trackTsize = dt;
            if isempty(tmp.trackTtickdist)
                tmp.trackTtickdist = dt/10;
            end         
        end
        
        function set.trackFsize(tmp,df)
            tmp.trackFsize = df;
            if isempty(tmp.trackFtickdist)
                tmp.trackFtickdist = 25e3;
            end            
        end
        
        function settracklim(tmp,n,m)
        if nargin <3
            m = n;
        end
        
        dt  = tmp.trackTsize;
        df  = tmp.trackFsize;        
        
        tmp.trackTlim       = [m-1 m]*dt;
        tmp.trackTtick      = (m-1)*dt:tmp.trackTtickdist:m*dt;
        tmp.trackTtickstr   = tmp.trackTtick;            
        
        if isempty(tmp.trackFlim)
            tmp.trackFlim       = [0 df];
        end
        tmp.trackFtick      = tmp.trackFlim(1):tmp.trackFtickdist:tmp.trackFlim(2);
        tmp.trackFtickstr   = tmp.trackFtick/1e3;
        end
        
        function settrack(tmp,n,m)
        if nargin ==3
            tmp.settracklim(n,m);
        end              
        h = tmp.trackhandle(n);
        set(h,'YLim',tmp.trackFlim,'YDir','normal');
        set(h,'XLim',tmp.trackTlim,'XDir','normal');
        set(h,'Ytick',tmp.trackFtick,'YTickLabel',tmp.trackFtickstr )
        set(h,'Xtick',tmp.trackTtick,'XTickLabel',tmp.trackTtickstr)
        set (h,'TickLength', [0.005 0])
            ylabel(h,'Freq. (kHz)');
            if n == tmp.trackNtotal
            xlabel (h,'Time (sec.)');            
            end
        
        tmp.trackNumber = n;        
        end
        
        
        function [x,y,I] = surf2image(tmp,X,Y,Z)
                h = tmp.trackhandle(tmp.trackNumber);
                set(h,'Units', 'pixels');
                pixpos = get(h,'Position');
                
                i = round(pixpos(3));
                j = round(pixpos(4));
                
                x  = linspace(tmp.trackTlim(1),tmp.trackTlim(2),i);
                ix = interval(x,[min(X) max(X)],'[]');
                
                y  = linspace(tmp.trackFlim(1),tmp.trackFlim(2),j);
                iy = interval(y,[min(Y) max(Y)],'[]');
                
                z = interpn(Y,X,Z,y(iy),x(ix)','linear');
                zdB =  pow2db(abs(z));                
                tmp.setClim(zdB);
                grayzdB = mat2gray(zdB,tmp.trackClim);
                [indzdB, ~] = gray2ind(grayzdB, length(tmp.trackCmap));
                I1 = ind2rgb(indzdB,tmp.trackCmap);
                
                I = ind2rgb(zeros(j,i), tmp.trackCmap);
                I(iy,ix,:) = I1;
            
        end
        
        function setClim(tmp,Z)
            switch tmp.trackClimMode
                case 'min'
                     z = Z(inf>Z & Z>-Inf);
                     tmp.trackClim = [min(z(:)) max(z(:))];
                case 'mean'
                    z = Z(inf>Z & Z>-Inf);
                     tmp.trackClim = [mean(z(:)) max(z(:))];
                otherwise
                    tmp.trackClim = [Z(1) Z(2)];
            end
        end
        
        function plottrack(tmp,P,F,T,n,m)
                tmp.trackNumber = n;
                if nargin <6
            	 m = n;
                end
                tmp.settracklim(n,m)
                h = tmp.trackhandle(n);
                [t,f,I] = surf2image(tmp,T,F,P);                
            % Image plot and set of axes
                axes(h)
                image(t,f,I);
                set (h,'Color',tmp.trackCmap(1,:,:))
                tmp.settrack(n)
        end
        
        function plotbox(tmp,X,Y,c,n,m)
                if nargin == 3
                    c =1;
                elseif nargin == 4
                    n = tmp.trackNumber;
                elseif nargin==5
                    m = n;
                    tmp.settracklim(n,m)
                elseif nargin==6
                    tmp.settracklim(n,m)
                end
                
                I = size(X,1);
                if ischar (c)
                    cbox = repmat(c,I,1);
                else
                    cbox = tmp.boxcolor(c,:);
                end
                
                h = tmp.trackhandle(n);                
                axes(h)
                for i = 1:I
                    x = X(i,:);
                    y = Y(i,:);
                    if tmp.trackTlim(1)<x(2) && x(1) < tmp.trackTlim(2)
                        if tmp.trackFlim(1)<y(2) && y(1) < tmp.trackFlim(2)
                        rectangle('Position',[x(1) y(1) diff(x) diff(y)],...
                                 'LineWidth',0.9,'LineStyle','-','EdgeColor',cbox(i,:));
                        
                        end
                    end
                end                
        end
        
        function plottext(tmp,x,y,tx,c,n,m)
                if      nargin == 4
                    c =1;
                elseif  nargin == 5
                    n = tmp.trackNumber;
                elseif  nargin == 6
                    m = n;
                    tmp.settracklim(n,m)
                elseif  nargin == 7
                    tmp.settracklim(n,m)
                end
                
                I = size(x,1);
                if ischar (c)
                    cbox = repmat(c,I,1);
                else
                    cbox = tmp.boxcolor(c,:);
                end
                h = tmp.trackhandle(n);                
                axes(h)
                for i = 1:I                    
                    if tmp.trackTlim(1)< x(i) && x(i) < tmp.trackTlim(2)
                        if tmp.trackFlim(1)<y(i) && y(i) < tmp.trackFlim(2)
                            text(x(i),y(i),tx(i,:),'Color',cbox(i,:),...
                            'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
                        end
                    end
                end
        end
        
        function plotpoints(tmp,x,y,c,n,m)
                if      nargin == 3
                    c = tmp.boxcolor;
                elseif  nargin == 4
                    n = tmp.trackNumber;
                elseif  nargin == 5
                    m = n;
                    tmp.settracklim(n,m)
                elseif  nargin == 6
                    tmp.settracklim(n,m)
                end
                
                if ischar (c)
                    if ~isletter(c)
                        cmarker = [ c tmp.boxcolor];
                    else
                        cmarker = [ '.' c];
                    end
                else
                    cmarker = '.';
                end
                
                I = size(x,1);
                h = tmp.trackhandle(n);                
                axes(h)
                hold on
                for i = 1:I                    
                    if tmp.trackTlim(1)< x(i) && x(i) < tmp.trackTlim(2)
                        if tmp.trackFlim(1)<y(i) && y(i) < tmp.trackFlim(2)
                            if ischar(c)
                                plot(h,x,y,cmarker);
                            else
                                hp = plot(h,x,y,cmarker);
                                set(hp,'MarkerEdgeColor',c);
                            end
                        end
                    end
                end
                hold off
        end
        
        function set.boxcolor(tmp,c)
            if ~ischar(c)
                tmp.boxcolor = lines(c);
            else
                tmp.boxcolor = c;
            end
        end
        
        function set.pageNote1(obj,sp)
            obj.pageNote1 = strcat('Fs = ',num2str(sp.fs/1e3,'%2.1f kHz'),...
                                                  ', wt = ',num2str(sp.wt*1e3,'%2.1f ms'),...
                                                  ', ov = ',num2str(sp.ov/sp.wn*100,'%%%2.2f'),...
                                                  ', df = ',num2str(sp.df,'%2.2f Hz'),...                                              
                                                  ', dt = ',num2str(sp.dt*1e3,'%2.3f ms'));
                                                  
        end
        
        function writenotes(tmp,h)
            
            if nargin <2
            h = tmp.pagehandle;
            end
            
            if ~verLessThan('matlab','8.3')
                h = figure(h);
            end
            %tmp.writetitle(h,1)           
                  
            pos = [tmp.framePos(1) 0.25 tmp.frameL 0.5];
            annotation(h,'textbox','Units','centimeters','position',pos,...
                          'String',tmp.pageNote1,...
                          'HorizontalAlignment','left','VerticalAlignment','bottom',...
                          'EdgeColor','none','FontSize',7.5)
            pos = [tmp.framePos(1) 0.23 tmp.frameL 0.5];         
            annotation(h,'textbox','Units','centimeters','position',pos,...
                          'String',tmp.pageNote2,...
                          'HorizontalAlignment','left','VerticalAlignment','top',...
                          'EdgeColor','none','FontSize',7.5)
            pos = [tmp.framePos(1)+tmp.framePos(3) tmp.framePos(2)+tmp.framePos(4) tmp.framel tmp.frameh];
            annotation(h,'textbox','Units','centimeters','position',pos,...
              'String',strcat(num2str(tmp.audiofileNumber),'/',num2str(tmp.audiofileNtotal)),...
              'HorizontalAlignment','center','VerticalAlignment','middle',...
              'EdgeColor','none','FontSize',14)
          
            pos = [tmp.framePos(1)+tmp.framePos(3) 0 tmp.framel tmp.frameh];
            annotation(h,'textbox','Units','centimeters','position',pos,...
                  'String',num2str(tmp.pageNumber),...
                  'HorizontalAlignment','center','VerticalAlignment','middle',...
                  'EdgeColor','none','FontSize',12,'FontWeight','normal')
        end
        
        function writetitle(tmp,h,n)
            
            if nargin ==2
                n   = tmp.trackNumber;                 
            end
            if nargin ==1
                h = tmp.pagehandle;
            end
            if isempty(tmp.trackTitleColor)
                c = 'k';
            else
                c = tmp.trackTitleColor;
            end
            tmp.pageTitle =  strcat('File ',num2str(tmp.audiofileNumber),': ',tmp.audiofileName);
            %pos = [tmp.framePos(1) tmp.framePos(1)+tmp.frameH tmp.frameL tmp.frameh];
            pos = tmp.trackPos(n,:);
            
            if ~verLessThan('matlab','8.3')
                h = figure(h);
            end
            
            annotation(h,'textbox','Units','centimeters','position',pos,...
                      'String',tmp.pageTitle,'Color',c,...
                      'HorizontalAlignment','left','VerticalAlignment','top',...
                      'EdgeColor','none','FontSize',14,'FontWeight','bold',...
                      'BackgroundColor','none','Interpreter','none')
        end

        function printpdf(tmp,h)
            
            if nargin <2
                h = tmp.pagehandle;
            end 
            drawnow
            if tmp. savekey
            set(h,'Color','w','PaperOrientation','Portrait','PaperUnits', 'centimeters',...
                              'PaperType','A4','PaperSize',tmp.pageSize,...
                              'PaperPosition', [0 0 tmp.pageSize])
            set(h,'PaperPositionMode', 'manual','InvertHardcopy','off')
            %saveas(h,tmp.savename,'pdf')
            %saveas(h,tmp.savename,'png')
            %print(h,tmp.savename,'-dpsc','-append','-loose','-r150','-painters')
            %print(h,tmp.savename,'-dpdf','-r150','-opengl')
            print(h,tmp.savename,'-dpdf','-r150','-painters')
            end
            end
        
        function newpage(tmp,h)
            if nargin <2
                h = tmp.pagehandle;
            end 
            clf(h);            
            tmp.newtracks(h);
            tmp.pageNumber = tmp.pageNumber+1; 
        end
        
        function clearpage(tmp,h)
            if nargin <2
                h = tmp.pagehandle;
            end 
            clf(h);         
        end
    end   
        
end

