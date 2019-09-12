function expression = getExpression(node,pattern,kineticMatrix,expression)
%--------------------------------------------------------------------------
% Get the expression given an enzyme form, patterns, and kineticMatrix
% 
% Inputs:        (node)  initial node index
%             (pattern)  specific pattern
%       (kineticMatrix)  stoichmetric matrix of the system, element are
%                        pseudo-first rate constants
%          (expression)  equation expression
%
% Outputs: (expression)  equation expression
%------------Pedro Saa 2014, adapted from Qi et al. 2009-------------------
%% 1. Use the singleNode to process the first node
[expression,pattern,nodes] = singleNode(node,pattern,kineticMatrix,expression);
node_new = nodes;
while (~isempty(pattern))
    nodes_old = node_new;
    node_new = [];
    for i = 1:length(nodes_old)
        node = nodes_old(i);
        [expression,pattern,nodes] = singleNode(node,pattern,kineticMatrix,expression);
        if(isempty(pattern))
            break;
        end
        node_new = [node_new,nodes];
    end
end

%% 2. Single node processing function
function [expression,pattern,nodes] = singleNode(node,pattern,kineticMatrix,expression)
nodes = []; expression_temp = []; number = 0; markLine = [];
for i = 1:size(pattern,1)        
    if (node == pattern(i,1))
        nodes = [nodes,pattern(i,2)];               
        number = number + 1;  
        markLine(number) = i;                      
        if(kineticMatrix{pattern(i,2),node}==0)
            expression_temp = [expression_temp; '%'];
        else
            expression_temp = [expression_temp;kineticMatrix(pattern(i,2),node)];
        end
    elseif (node == pattern(i,2))
        nodes = [nodes, pattern(i,1)];             
        number = number + 1;
        markLine(number) = i;                     
        if(kineticMatrix{pattern(i,1),node}==0)
            expression_temp = [expression_temp; '%'];
        else
            expression_temp = [expression_temp;kineticMatrix(pattern(i,1),node)];
        end
    end    
end
% check if get the new expression term
if(~isempty(expression_temp))    
    if (isempty(expression))       
        if (strfind(char(expression_temp(1)), '+'))
            expression = [expression '(' char(expression_temp(1)) ')' ];
        else
            expression = [expression char(expression_temp(1))];
        end
    else               
        if (strfind(char(expression_temp(1)), '+')) 
            expression = [expression '*(' char(expression_temp(1)) ')' ];
        else
            expression = [expression '*' char(expression_temp(1))];
        end
    end
    % if has more than one term
    if length(expression_temp) > 1
        for i = 2 : length(expression_temp)                       
            if (strfind(char(expression_temp(i)), '+'))
                expression = [expression '*(' char(expression_temp(i)) ')' ];
            else
                expression = [expression '*' char(expression_temp(i))];
            end
        end
    end
end
markLine = sort(markLine,'descend');
for i = 1:length(markLine)
    pattern(markLine(i),:) = [];
end