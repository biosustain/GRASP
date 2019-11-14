function [nodeList,edge,kineticMatrix,forwardFlux] = readInput(filename)
% Reads input file of the enzymatic reaction pattern.
% 
% Based on https://doi.org/10.1186/1471-2105-10-238 
%
%
% USAGE:
%
%    [nodeList, edge, kineticMatrix, forwardFlux] = readInput(filename)
%
% INPUT:
%    filename (char):   pattern name
% 
% OUTPUT:
%    nodeList (int vector):      list with the number of nodes
%    edge (int cell):            edge list
%    kineticMatrix (cell):       stoichmetric matrix of the system, element are pseudo-first rate constants
%    forwardFlux (int matrix):   list with the edges of the forward reactions
%
% .. Authors:
%       - Pedro Saa     2016 original code, adapted from Qi et al. 2009

% 1. Retrieve information from input file
currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
filepath = fullfile(currentPath{1}, '..', '..', 'patterns', filename);

try
    InputFileFlow = textread([filepath,'.txt'],'%s');
catch
    error(['The file for the mechanism named ', filename, ' could not be found: ', filepath, ...
            newline, ...
           'Please make sure the kinetic mechanism name corresponds to a file ', ...
           'with the same name in the folder patterns. Also the file extension should be ".txt"']);
end

FileLength = length(InputFileFlow); nodeList = [];
for i = 1 : FileLength/3     
    num1 = str2num(InputFileFlow{i*3-2});
    num2 = str2num(InputFileFlow{i*3-1});
    if(isempty(num1) || isempty(num2))
        errordlg(['No right information in input file!';' Please select another input file! '],'File Open Error');
        edge = {};
        kineticMatrix = {};
        return;
    else
    nodeList = [nodeList,num1];
    nodeList = [nodeList,num2];
    end
end

% 2. Recover forward elementary fluxes
forwardFlux = [];
for k = 1:length(nodeList)/4
    forwardFlux = [forwardFlux;nodeList(4*k-3),nodeList(4*k-2)];
end

% 3. Build the kineticMatrix from input file
nodeList = unique(nodeList);
kineticMatrix{max(nodeList),max(nodeList)} = [];
for i = 1:FileLength/3     
    num1 = str2num(InputFileFlow{i*3-2});
    num2 = str2num(InputFileFlow{i*3-1});
    if(isempty(kineticMatrix{num1,num2}))
        kineticMatrix{num1, num2} = InputFileFlow{i*3};
    else
        kineticMatrix{num1, num2} = [kineticMatrix{num1, num2},'+',InputFileFlow{i*3}]; 
    end
end

% 4. Set up the edge information
edgeMatrix = zeros(length(nodeList),length(nodeList));
edge       = cell(size(kineticMatrix,1),3);
edgeIndex  = 0;
for i = 1:size(kineticMatrix,1)
    edgeNumber = 0;
    edgeList   = [];
    for j = 1 : size(kineticMatrix,1)
        if(i == j || (isempty(kineticMatrix{i,j}) && isempty(kineticMatrix{j,i})))  
            continue;
        else
            edgeNumber = edgeNumber + 1;
            edgeList   = [edgeList,j];
            if(j>i) 
                edgeIndex = edgeIndex + 1;
                edgeMatrix(i,j) = edgeIndex;
                edgeMatrix(j,i) = edgeIndex;
            end
        end
    end
    edge(i,:) = {edgeNumber,edgeList,edgeMatrix(i,edgeList(:))};
end

% 5. Set [] to 0 in the final kineticMatrix
for i = 1:size(kineticMatrix,1)
    for j = 1:size(kineticMatrix,1)
        if(isempty(kineticMatrix{i,j}))
            kineticMatrix{i,j} = 0;
        end
    end
end