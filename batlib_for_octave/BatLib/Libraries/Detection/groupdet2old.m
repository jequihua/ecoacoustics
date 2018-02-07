function [G,dX] = groupdet2(X,m,gmin,gmax,dmax)

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
    m = 3;
end

if ~issorted(X)
    disp('Warning: detections are not sorted')
end

if m<3
    m=3;
end

  I    = length(X);
  G    = zeros(I,1);
  dX   = zeros(I,1);
  
    if I >2

     D = squareform(pdist(X,'euclidean'));
     D(D>dmax) = Inf;

      for i = 1:I
          if G(i)== 0          
              g = max(G)+1;
              S = 1;
              j = i+1;
              while j <=I && S
                  if G(j)==0
                    d1  = D(i,j);
                    k = j+1;
                    while k <=I && S 
                        if G(k)==0
                        d2 = D(j,k);
                        d  = d2/d1;
                            if gmin <= d && d <=gmax
                            G(i) = g;
                            G(j) = g;
                            G(k) = g;
                            S    = 0;
                            dX(i) = d1;
                            dX(j) = d1;
                            dX(k) = d2;
                            end
                        end
                        k = k+1;
                    end
                  end
                  j=j+1;
              end

              if S
                  G(i) = g;              
              end
          else
              g     = G(i);
              ig    = find(G==g);
              if i == ig(end)
                  d1 = dmean(D,ig,m);
                  S    = 1;
                  j   = i+1;
                  while j <=I && S
                      d2 = D(i,j);
                      d  = d2/d1;
                      if gmin <= d && d <=gmax                      
                          G(j) = g;
                          dX(i) = d1;
                          dX(j) = d2;
                          S    = 0;
                      end
                      j=j+1;
                  end
              end
          end 
      end
    else
        G = [1;2];
    end
    
 
    function d = dmean(D,ind,m)        
         n     = length(ind);          
         dbag  = nan(m-1,1);
         for i = 1:min(n-1,m-1)
                j   = n-i+1;
                dbag(i,1)  = D(ind(j-1),ind(j));
         end
         d     = nanmean(dbag); 
  

