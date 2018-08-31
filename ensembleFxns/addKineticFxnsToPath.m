function addKineticFxnsToPath(ensemble)
% Add kinetic fxns to the path
% 
% ------------------ Pedro Saa, Marta Matos 2018 ---------------------------------------
%
rootName = 'reactions';
for ix = 1:ensemble.numStruct
    tmpName = [rootName,num2str(ix)];
    addpath(tmpName);
end