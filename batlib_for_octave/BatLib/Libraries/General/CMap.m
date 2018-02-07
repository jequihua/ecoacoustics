classdef CMap < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N0     = 16;
        N      = 128;
        cmap   = 'jet';
        cfun   = 'linear';
        cvalue = 1;
        gamma = 2;
        tol   = 0.5;
        char0 = [0 17 34 51 68 85 102 119 136 153 170 187 204 221 238 255]/255;
        char1 = [2 3 4 6 8 11 14 19 25 33 42 55 72 100 150 255]/255;
        char2 = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]/255;
        char3 = [3 7 13 18 23 30 38 48 65 87 120 150 190 230 250 255]/255;
        char4 = [7 10 15 20 28 40 54 75 93 120 145 175 205 230 254 255]/255;
        char5 = [20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95]/255;
        
        cbri =      [0     0     0
                     7     0   145
                     0    46   202
                     0   120   220
                     0   160   255
                     0   210   255
                     2   220   205
                    30   230   170
                    70   240   110
                   100   255    84
                   223   255     0
                   255   245     0
                   255   220     0
                   255   160     0
                   250    80     0
                   230     0     0]/255;
               
        cdar =      [0     0     0
                     7     0   145
                     0    46   202
                     0   120   220
                     0   160   255
                     0   210   255
                     2   220   205
                    30   230   170
                    70   240   110
                   100   255    84
                   223   255     0
                   255   245     0
                   255   220     0
                   255   160     0
                   250    80     0
                   230     0     0]/255;     
    end
    
    methods
        function c = colorbri(tmp,varargin)
            [N,cfun,cvalue,~] = tmp.getvariables(varargin);
            c0 = tmp.cbri;
            c  = tmp.Cfun(c0,N,cfun,cvalue);            
        end        
        function c = colordar(tmp,varargin)
            [N,cfun,cvalue,~] = tmp.getvariables(varargin);
            c0 = tmp.cdar;
            c = tmp.Cfun(c0,N,cfun,cvalue);           
        end        
        function c = colormap(tmp,varargin)
            [N,cfun,cvalue,cmap] = tmp.getvariables(varargin);                       
            if isstr(cmap)
                try
                    eval(['c0=' cmap '(' num2str(N) ')']);
                end
                try
                    c0 = tmp.(cmap)(N,cfun,cvalue);
                end
            elseif size(cmap,2)==3;
                    c0  = cmap;                    
            end            
            c = tmp.Cfun(c0,N,cfun,cvalue);
        end        
        function c = grayneg(tmp,varargin)
            [N,cfun,cvalue,~] = tmp.getvariables(varargin);
            c0 = 1-gray(N);
            c = tmp.Cfun(c0,N,cfun,cvalue);
        end
        function c = gray(tmp,varargin)
            [N,cfun,cvalue,~] = tmp.getvariables(varargin);
            c0 = gray(N);
            c = tmp.Cfun(c0,N,cfun,cvalue);
        end        
        function c = jet(tmp,varargin)
            [N,cfun,cvalue,~] = tmp.getvariables(varargin);
            c0 = jet(N);
            c = tmp.Cfun(c0,N,cfun,cvalue);
        end
        function [N,cfun,cvalue,cmap] = getvariables(tmp,x)
            if length(x) == 0
                N      = tmp.N;
                C      = tmp.cfun;
                cvalue = tmp.cvalue;
                cmap   = tmp.cmap; 
            elseif length(x) == 1
                if isscalar(x{1})
                    N      = x{1}; 
                    C      = tmp.cfun;
                    cvalue = tmp.cvalue;
                    cmap   = tmp.cmap;
                elseif  size(x,2)== 3;
                    N      = size(x{1},1); 
                    C      = tmp.cfun;
                    cvalue = tmp.cvalue;
                    cmap   = x{1};
                end
            elseif length(x) == 2
                if isscalar(x{1})
                    N      = x{1}; 
                    C      = x{2};
                    cvalue = tmp.cvalue;
                    cmap   = tmp.cmap;
                elseif  size(x,2)== 3;
                    N      = size(x{1},1);
                    C      = x{2};
                    cvalue = tmp.cvalue;
                    cmap   = x{1};                    
                end
            elseif length(x) == 3
                cmap   = tmp.cmap;
                N      = x{1};
                C      = x{2};
                cvalue = x{3};                                
            elseif length(x) == 4
                cmap   = x{1};
                N      = x{2};
                C      = x{3};
                cvalue = x{4};                 
            end
            
            
            if ischar(C)
                cfun   = C;
            elseif isvector(C)
                cvalue = C;
                cfun  = 'gamma';     
            end
        end
        function c = Cfun(tmp,c0,N,cfun,cvalue)            
            if size(c0,1)~=N;
                N0 = size(c0,1);
                [y,y0] = tmp.getCfun(N0,N,cfun,cvalue);
                c  = zeros(N,3);            
                for i = 1:3, 
                    c(:,i) = interpn(y0,c0(:,i),y,'pchip'); 
                end
            else
                switch cfun
                    case 'gamma'
                        c = imadjust(c0,[0 1],[],cvalue);
                    case 'linear'
                        c = c0;  
                    case {'char0','char1','char2','char3','char4','char5'}
                        [y,y0] = tmp.getCfun(N,N,cfun,cvalue);
                        c  = zeros(N,3);            
                        for i = 1:3, 
                            c(:,i) = interpn(y0,c0(:,i),y,'pchip'); 
                        end     
                    end
                
            end
        end            
        function [y,y0] = getCfun(tmp,N0,N,cfun,cvalue)          
            switch cfun
                    case 'gamma'
                        x0 = linspace(0,1,N0);
                        y0 = imadjust([x0' x0' x0'],[0 1],[],cvalue);
                        y0 = y0(:,1);
                        x = linspace(0,1,N); 
                        y = imadjust([x' x' x'],[0 1],[],cvalue);
                        y = y(:,1);
                    case 'linear'
                        y0 = linspace(0,1,N0);
                        y  = linspace(0,1,N);                        
                    case {'char0','char1','char2','char3','char4','char5'}
                       x0 = tmp.(cfun);
                       y0 = linspace(0,1,tmp.N0);
                       
                       x = linspace(0,1,N);
                       y = pchip(x0,y0,x);
                       y = max(0,min(y,1));
                       if  N == N0;
                           y0 = x;
                       end
                end
        end
    end
end

