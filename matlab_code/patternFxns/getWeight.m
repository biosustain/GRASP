function W = getWeight(shortPath,forwardFlux)
% Function used re-weight a path in order to obtain alternative paths.
%
%
% USAGE:
%
%    W = getWeight(shortPath, forwardFlux)
%
% INPUT:
%    shortPath (int vector):      common path
%    forwardFlux (int matrix):	  list with the edges of the forward reactions
% 
% OUTPUT:
%    W (int vector):	weight vector
%
% .. Authors:
%       - Pedro Saa     2014 original code

%% 1. Define edges to re-weight
idx = []; pathDist = length(shortPath)-1;
for i = 1:size(forwardFlux,1)
    for u = 1:pathDist
        if  isequal(forwardFlux(i,:),[shortPath(u),shortPath(u+1)]) == 1
            idx = [idx;i];
        end
    end
end
%% 2. Assing new weights to the edges
W = ones(size(forwardFlux,1),1);
W(idx) = 2; % 2-times the usual weight (this works in MATLAB 2013)
% For previous MATLAB versions, used W(ix) = 0.5, instead