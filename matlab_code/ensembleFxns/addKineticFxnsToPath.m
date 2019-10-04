function addKineticFxnsToPath(ensemble)
% Adds all the model's kinetic functions to Matlab's path.
%
%
% USAGE:
%
%    addKineticFxnsToPath(ensemble)
%
% INPUT:
%    ensemble (struct):   model ensemble
%
% .. Authors:
%       - Pedro Saa      2016 original code
%       - Marta Matos    2018 rootName based on current directory

currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
rootName =  fullfile(currentPath{1}, '..', '..', 'reactions', strcat(ensemble.description, '_'));

for ix = 1:ensemble.numStruct
    tmpName = [rootName,num2str(ix)];
    addpath(tmpName);
end