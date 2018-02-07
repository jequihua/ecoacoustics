function filepath = changepath(filepath,varargin)

if length(varargin)==1
    path      = varargin{1};
    system1   = path.Sytem;
    system2   = computer; 
    WinHeader = path.WinHeader;
    LinHeader = path.LinHeader;
    
elseif length(varargin)== 2
    path      = varargin{1};
    system1   = varargin{2};
    system2   = path.Sytem;        
    WinHeader = path.WinHeader;
    LinHeader = path.LinHeader;
elseif length(varargin)== 3
    system1    = computer;
    system2    = varargin{1};
    WinHeader =  varargin{3};
    LinHeader =  varargin{4};
elseif length(varargin)==4
    system1    = varargin{1};
    system2    = varargin{2};
    WinHeader =  varargin{3};
    LinHeader =  varargin{4};
end
        
switch system2
        case {'PCWIN','PCWIN64'}
                switch system1
                 case {'GLNXA64'}
                        filepath = strrep(filepath,LinHeader,WinHeader);

                        oldSubstr = '/';
                        newSubstr = '\';
                        filepath = strrep(filepath, oldSubstr, newSubstr);
                end

            
        case {'GLNXA64'}
                    switch system1
                    case {'PCWIN','PCWIN64'}
                        filepath = strrep(filepath,WinHeader,LinHeader);

                        oldSubstr = '\';
                        newSubstr = '/';
                        filepath = strrep(filepath, oldSubstr, newSubstr);
                    end
end