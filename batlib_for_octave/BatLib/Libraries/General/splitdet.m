function [M,iM] = splitdet(Y)
  
  N = nzmax(Y);
  M  = zeros(N,1);
  iM = zeros(N,2);
  y = zeros(N+1,1);
  y(1:N) = Y;
  
  m = 0;
  g = 0;
  z = 0;
  
  for i = 1:N+1;
      if     y(i)==1 && z==0
          z = 1;
          g = g+1;
          m = m+1;
          iM(g,1)=i;
      elseif y(i)==0 && z==1
          M(g)= m;
          iM(g,2)= i-1;
          z = 0;
          m = 0;
      elseif y(i)==1 && z==1
          m = m+1;
      elseif y(i)==0 && z==0                              
      end                          
  end
  M = M(1:g);
 iM = iM(1:g,1:2);