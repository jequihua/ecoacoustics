function [G,dX] = groupdet2(X,m,gmin,gmax,dmax)

if nargin<5
    dmax(1) = 250e-3;
    dmax(2) = 10e3;
end
if nargin<4
    gmax(1) = 1.5;
    gmax(2) = 1.15;
end
if nargin<3
    gmin(1) = 0.5;
    gmin(2) = 0.85;
end
if nargin<2
    m = 2;
end

if ~issorted(X(:,1))
    disp('Warning: detections are not sorted')
end

if m<2
    m=2;
end

  I    = size(X,1);
  G    = zeros(I,1);
  dX   = zeros(I,2);
  
    if I >2
     while any (~G)         
         i = find(~G,1,'first');
         g = max(G)+1;
         key = 1;
         while key              
              key1 = 0;
              for j= ifilt2(X,i,I,gmin,gmax,dmax)
                  if G(j)==0
                    d1  = dist(X,i,j,dmax);                    
                    for k = ifilt2(X,j,I,gmin,gmax,dmax)
                        if G(k)==0
                        d2 = dist(X,j,k,dmax);
                        d  = d2/d1;
                            if gmin(1) <= d && d <=gmax(1)
                            G(i) = g;
                            G(j) = g;
                            G(k) = g;
                            dX(i) = d1;
                            dX(j) = d1;
                            dX(k) = d2;
                            key1  = 1;
                            break
                            end
                        end                        
                    end
                    if key1
                        break
                    end
                  end
              end
              
           if key1               
               while j<I             
                   ig  = find(G==g)';
                   j   = ig(end);                   
                   d1  = dmean(X,ig,m,dmax);
                   key2 = 1; 
                      for k = ifilt2(X,j,I,gmin,gmax,dmax)  
                          if G(k)==0
                            d2 = dist(X,j,k,dmax);
                            d  = d2/d1;
                            if gmin(1) <= d && d <=gmax(1)
                              G(j) = g;
                              G(k) = g;
                              dX(j) = d1;
                              dX(k) = d2;
                              key2  = 0;
                              break
                            end
                          end                                
                      end
                      if key2
                          key  = 0;
                          break
                      end                    
               end
           end           
           if key
           G(i) = g;
           key  = 0;
           end              
         end
     end
    else
        G = [1;2];
    end
    
 
    function d = dmean(X,ind,m,dmax)        
         n     = length(ind);          
         dbag  = nan(m-1,1);
         for i = 1:min(n-1,m-1)
                j   = n-i+1;
                dbag(i,1)  = dist(X,ind(j-1),ind(j),dmax);                
         end
         d     = nanmean(dbag,1);
         
        function d = dist(X,i,j,dmax)
            d = X(j,1)-X(i,1);            
            if d > dmax(1)
                d = Inf;
            end
            
            
            function j = ifilt2(X,i,I,gmin,gmax,dmax)
                j = false(1,I);
                    for k = i+1:I
                        x  = X(k,2)/X(i,2);
                        dX = abs(X(k,2)-X(i,2)); 
                        if dX <= dmax(2)
                            if gmin(2) <= x && x <=gmax(2)
                                j(k)=1;
                            end
                        end
                    end
                j = find(j);

                
                
              

