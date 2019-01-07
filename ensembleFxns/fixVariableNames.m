function outputCell = fixVariableNames(cellWithVariables, prefix)
%--------------------------------------------------------------------------
% Takes in a cell where the first column are metabolite or reaction names
% and adds the specified prefix to each name.
%
% Inputs:       cellWithVariables (cell), prefix (string)
%
% Outputs:      outputCell (cell)
%--------------------- Marta Matos 2019 -----------------------------------


for row=2:size(cellWithVariables,1)
    cellWithVariables{row,1} = strcat(prefix, '_', cellWithVariables{row,1});
end

end

