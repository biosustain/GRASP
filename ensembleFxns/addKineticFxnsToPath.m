function addKineticFxnsToPath(ensemble)
% Add kinetic fxns to the path
% 
% ------------------ Pedro Saa, Marta Matos 2018 ---------------------------------------
%
rootName =  strcat('reactions_', ensemble.description, '_');

for ix = 1:ensemble.numStruct
    tmpName = [rootName,num2str(ix)];
    addpath(tmpName);
end