function plotcall(P,F,T,S,Pi,Fi,Ti,i0,j0,e,f,t,stft,m)
        
        cmap = jet;

        h1 = figure(1);
        clf(h1)
        title(m)
        hold on
            surf(Ti,Fi,10*log10(Pi),'edgecolor','none'); axis tight; view(0,90);
            colormap(cmap)
            plot3(T(j0),F(i0),10*log10(P(i0,j0)),'*w');

            zlim([-200 0]) 
            x = stft.tplim(m,1);
            y = stft.fplim(m,1);
            dx = diff(stft.tplim(m,:));
            dy = diff(stft.fplim(m,:));
            rectangle('position',[x y dx dy],'EdgeColor','r', 'LineWidth',2)

            x = T(1);
            y = F(1);
            dx = T(end)-T(1);
            dy = F(end)-F(1);
             rectangle('position',[x y dx dy],'EdgeColor','w', 'LineWidth',2)
         hold off
         
        figure(2);
        clf(2)
        
        subplot(2,1,1)
            hold on
            surf(T,F,10*log10(P.*S),'edgecolor','none'); 
            plot3(t,f,10*log10(e),'-k','LineWidth',2);
            hs = get(h1,'Children');
            caxis(get(hs,'CLim'));
            zlim(gca,get(hs,'ZLim'));
            xlim(gca,get(hs,'XLim'));
            ylim(gca,get(hs,'YLim'));
            set(gca,'Color',cmap(1,:,:))        
             xlabel(nnz(S)/nzmax(S)*stft.aframe)
             axis tight; view(0,90);
             hold off
         subplot(2,1,2)
             df = diff(f,1,2);
             dt = diff(t,1,2);
             plot(t(1:end-1),df./dt,'-k','LineWidth',2);
             xlim(gca,get(hs,'XLim'));
             ylim(gca,[-1 1]*2e7);
            %set(gca,'Color',cmap(1,:,:))    
            if nnz(S)/nzmax(S)>stft.amax
                pause(1)
            else
                pause(0.05)
            end
            

                
         
         
        %plot3(T,F,E,'-w');
        