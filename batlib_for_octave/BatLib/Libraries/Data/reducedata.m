function Obj = reducedata(Obj,iN)

    VarNames = fieldnames(Obj);
    
    for i = 1: length(VarNames)        
        Variable = VarNames{i};
        [N,M] = size(Obj.(Variable));
        if N == length(iN)
            if islogical(Obj.(Variable))
                
                Obj.(Variable)(iN,:) = false(nnz(iN),M);
                
            elseif isnumeric(Obj.(Variable))
                
                %Obj.(Variable)(iN,:) = nan(nnz(iN),M);
                
            elseif iscell(Obj.(Variable))
                
                Obj.(Variable)(iN,:) = cell(nnz(iN),M);
                
            end
        end
        
        if isstruct(Obj.(Variable))
            disp(Variable)
             Obj.(Variable) = reducedata(Obj.(Variable),iN);
        end
    end
end

