
function [G,dX] = groupdet1(X,m,gmin,gmax,dmax)

if nargin<5
    dmax = 250e-3;
end
if nargin<4
    gmax = 1.5; 
end
if nargin<3
    gmin = 0.5; 
end
if nargin<2
    m = 2;
end

if ~issorted(X)
    disp('Warning: detections are not sorted')
end

if m<2
    m=2;
end

  I    = length(X);
  G    = zeros(I,1);
  dX   = zeros(I,1);
  
    if I >2
     while any (~G)         
         i = find(~G,1,'first');
         g = max(G)+1;
         key = 1;
         while key              
              key1 = 0;
              for j=i+1:I
                  if G(j)==0
                    d1  = dist(X,i,j,dmax);                    
                    for k = j+1:I
                        if G(k)==0
                        d2 = dist(X,j,k,dmax);
                        d  = d2/d1;
                            if gmin <= d && d <=gmax
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
                      for k = j+1:I  
                          if G(k)==0
                            d2 = dist(X,j,k,dmax);
                            d  = d2/d1;
                            if gmin <= d && d <=gmax
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
                dbag(i,:)  = dist(X,ind(j-1),ind(j),dmax);
         end
         d     = nanmean(dbag);
         
        function d = dist(X,i,j,dmax)
            d = abs(X(j)-X(i));
            if d > dmax
                d = Inf;
            end


