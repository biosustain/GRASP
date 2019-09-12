function addKineticFxnsToPath(ensemble)
% Add kinetic fxns to the path
% 
% ------------------ Pedro Saa, Marta Matos 2018 ---------------------------------------
%

currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
rootName =  fullfile(currentPath{1}, '..', '..', 'reactions', strcat(ensemble.description, '_'));

for ix = 1:ensemble.numStruct
    tmpName = [rootName,num2str(ix)];
    addpath(tmpName);
end