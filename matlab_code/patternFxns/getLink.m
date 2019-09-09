function [edgeMatrix,edgeList] = getLink(edge)
%--------------------------------------------------------------------------
% Get the edgeMatrix and edgeList data structures
% 
% Inputs:         (edge)  edge information structure
%
% Outputs:  (edgeMatrix)  matrix with the edge index as the elements
%             (edgeList)  edge list
%------------Pedro Saa 2014, adapted from Qi et al. 2009-------------------

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