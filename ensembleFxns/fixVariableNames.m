function outputCell = fixVariableNames(cellWithVariables, prefix)
%FIXVARIABLENAMES Summary of this function goes here
%   Detailed explanation goes here


disp(cellWithVariables);

for row=2:size(cellWithVariables,1)
    cellWithVariables{row,1} = strcat(prefix, '_', cellWithVariables{row,1});
end

disp(cellWithVariables);

end

