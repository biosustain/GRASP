function [revMatrix,forwardFlux,metList] = reactionPattern(patternName,reactionName,flag,strucIdx, promiscuousRxnI)
% Main function to build the reaction rate equation and to extract the
% required parameters for the subsequent sampling.
%
%
% USAGE:
%
%    [revMatrix, forwardFlux, metList] = reactionPattern(patternName, reactionName, flag, strucIdx, promiscuousRxnI)
%
% INPUT:
%    patternName (char vector):   input file which descrips the enzyme mechanism
%    reactionName (char):         name of the reaction
%    flag (int):                  if 1, write the reaction rate function (optional) else if 2, write the corresponding allosteric mechanism
%    strucIdx (int):              ID of the model structure
%    promiscuousRxnI (int):       index of promiscuous reaction
% 
% OUTPUT:
%    revMatrix (int vector):    reversibility matrix with all the cycles in the reaction pattern
%    forwardFlux (int matrix):  list with the edges of the forward reactions
%    metList (char cell):       list with the reaction metabolites
%
% .. Authors:
%       - Pedro Saa         2014 original code
%       - Marta Matos       2018 extended for promiscuous reactions
%       - Nicholas Cowie	2019 [TODO what were you doing?]

% 1. Read the input file and extract the required information
[Nodes,Link,KineticMatrix,forwardFlux] = readInput(patternName);
[LinkMatrix,LinkList] = getLink(Link);
prodIndex = [];
subsIndex = [];
prodNum = [];
counter = 1; 
for i = 1:size(KineticMatrix,1)
    
    for j = 1:size(KineticMatrix,2)
        
        % extract indices for reactions producing P, Q or R
        if (isempty(strfind(KineticMatrix{i,j},'P'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'Q'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'R'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'S'))~=1)
            
            prodIndex = [prodIndex;j,i];
            
            % indexes of P are saved for subsequent operations
            if  (isempty(strfind(KineticMatrix{i,j},'P'))~=1)
                prodNum = [prodNum;j,i];
                numTerm{counter,1} = ['+',KineticMatrix{j,i}];
                numTerm{counter,2} = ['-',KineticMatrix{i,j}];
                counter = counter + 1;
            end
            
        % extract indices for reactions consuming A, B or C
        elseif  (isempty(strfind(KineticMatrix{i,j},'A'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'B'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'C'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'D'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'E'))~=1 ||...
                isempty(strfind(KineticMatrix{i,j},'F'))~=1)
            subsIndex = [subsIndex;i,j];
        end
    end
end

tempPath  = colVector(cycledPaths(unique(subsIndex(:,1)),unique(prodIndex(:,1)),forwardFlux));
tempPath = tempPath(~cellfun('isempty',tempPath)); % get a column vector with all paths
revMatrix = zeros(size(tempPath,1),size(forwardFlux,1));

for k = 1:size(tempPath,1)
    tempPathList = [];
    
    % Extract nodes from the paths found
    for j = 1:size(tempPath{k},2)-1
        tempProd = isequal(tempPath{k}(j),tempPath{k}(j+1));
        if tempProd ~= 1
            tempPathList = [tempPathList;[tempPath{k}(j),tempPath{k}(j+1)]];
        else
            tempPathList = [tempPathList;zeros(1,2)];
        end
    end
    for u = 1:size(forwardFlux,1)
        for i = 1:size(tempPathList,1)
            
            % Compare nodes directionally
            if (isequal(forwardFlux(u,:),tempPathList(i,:))==1)
                revMatrix(k,u) = 1;
                break
            end
        end
    end
end

revMatrix = unique(revMatrix,'rows');

% 2. Process using alpha-numerical algebra
NodeNumber = length(Nodes);
tempPath   = [];
MaxIndex   = find(cell2mat(Link(:,1))==max(cell2mat(Link(:,1))));
for i = 1:NodeNumber
    if i == MaxIndex(1)
        continue;
    else
        tempPath = algebra(tempPath,Link{i,3}');
    end
end
Pattern = tempPath;       % all possible pattern list

% 3. Calculate every pattern expression corresponding to each state
PatternNumber = length(Pattern);
pathway = [];
for i = 1:PatternNumber
    try
        pattern = Pattern(i,:)';
    catch error
        error('There is  a good chance your mechanism has an error in the enzyme states numbering.');
    end
    
    pathway_tempPath = [];
    for j = 1:length(pattern)
        pathway_tempPath = [pathway_tempPath;LinkList(pattern(j),:)];
    end
    pathway(:,:,i) = pathway_tempPath;
end
expression = cell(length(Nodes),PatternNumber);
for i = 1:length(Nodes)
    node = Nodes(i);
    for j = 1:PatternNumber
        exp = '';
        exp = getExpression(node,pathway(:,:,j),KineticMatrix,exp);
        expression{i,j} = exp;
    end
end

% 4. Output fraction and denominator
state       = cell(length(Nodes),1);
Denominator = 'D = ';
for i = 1:length(Nodes);
    tempPath = char(expression{i,1});
    if (~isempty(strfind(tempPath,'%')))
        state{i} = [char(state{i}) ''];
    else
        state{i} = [char(state{i}) ' ' char(expression{i, 1})];
    end
    for j = 2 : PatternNumber
        tempPath = char(expression{i,j});
        if(~isempty(strfind(tempPath,'%')))
            state{i} = [char(state{i}) ''];
        else
            if(isempty(char(state{i})));
                state{i} = [char(state{i}) char(expression{i,j})];
            else
                state{i} = [char(state{i}) '+' char(expression{i,j})];
            end
        end
    end
end
Denominator = [Denominator char(state{1})]; state{1} = [char(state{1}) ';'];
for i = 2: length(Nodes);
    Denominator = [char(Denominator) '+' char(state{i})];
    state{i} = [char(state{i}) ';'];
end
if (length(strfind(Denominator,'+')) == length(Denominator)-3)
    Denominator = '';
end
Denominator = [Denominator ';'];
for i = 1 : length(Nodes);
    state{i} = ['E' num2str(i) '=' char(state{i})];
end

% 5. Give the rate constant list and the concentration list
k = 1;
for i = 1:length(Nodes)
    for j = 1:length(Nodes)
        if (LinkMatrix(i,j)~=0 && ~isempty(KineticMatrix{i,j}))
            tempPath = KineticMatrix{i, j};
            if(isnumeric(tempPath))
                rateList(k,1) = KineticMatrix(i,j);
                k = k+1;
            else
                RatetempPath = regexp(tempPath,'+','split');
                for m = 1:length(RatetempPath)
                    rateList(k,1) = RatetempPath(m);
                    k = k+1;
                end
            end
        end
    end
end
metList = {};
j = 1; k = 1;
for i = 1:length(rateList)
    if (rateList{i} ~= 0)
        tempPath = char(rateList{i});
        if(strfind(tempPath,'.*'))
            metList{k} = tempPath(strfind(tempPath,'.*')+1:end);
            k = k+1;
            rateList_new{j} = tempPath(1:strfind(tempPath,'.*')-1);
            j = j+1;
        else
            rateList_new{j} = tempPath;
            j = j+1;
        end
    end
end
rateList = sort(rateList_new);
rateList = unique(rateList);
metList = unique(metList);

% 6. Output the .m file (optional)
reactionName = [reactionName,num2str(strucIdx)];
if (flag == 1)                   % Non-allosteric reaction
    buildReaction(state,rateList,metList,numTerm,prodNum,reactionName, promiscuousRxnI);
elseif (flag == 2)               % Allosteric reaction
    buildReaction(state,rateList,metList,numTerm,prodNum,[reactionName,'Catalytic'], promiscuousRxnI);
end