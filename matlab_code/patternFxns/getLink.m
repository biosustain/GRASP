function [edgeMatrix,edgeList] = getLink(edge)
% Get the edgeMatrix and edgeList data structures.
%
% Based on https://doi.org/10.1186/1471-2105-10-238 
%
%
% USAGE:
%
%    [edgeMatrix, edgeList] = getLink(edge)
%
% INPUT:
%    edge (int cell):     edge information structure
% 
% OUTPUT:
%    edgeMatrix (int matrix):  matrix with the edge index as the elements
%    edgeList (int matrix):    edge list
%
% .. Authors:
%       - Pedro Saa     2014 original code, adapted from Qi et al. 2009

%% 1. Read edge information
nodeNum = size(edge,1);
edgeMatrix = zeros(nodeNum,nodeNum);

%% 2. Build edge matrix and edge list 
for i = 1:nodeNum
    neighborNumber = length(edge{i,2});
    for j = 1:neighborNumber
        edgeMatrix(i,edge{i,2}(j)) = edge{i,3}(j);
        edgeMatrix(edge{i,2}(j),i) = edge{i,3}(j);
    end
end
edgeList = [];
for i  = 1:nodeNum
    for j = i:nodeNum
        if (edgeMatrix(i,j)~=0)
            edge = sort([i j]);
            edgeList = [edgeList;edge];
        end
    end
end