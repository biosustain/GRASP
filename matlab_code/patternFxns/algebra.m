function A = algebra(B,C)
% Calculate wang algebra result for two sets of patterns
%
% Based on https://doi.org/10.1186/1471-2105-10-238 
%
% USAGE:
%
%    A = algebra(B, C)
%
% INPUTS:
%    B (`struct`):	input combination of link-number list string cell
%    C (`double`):  input of number string cell (only one link)
%
% OUTPUT:
%    A (`struct`):  Wang Algebra result calculate from input B and C, 
%                   result is a string cell
%
% .. Authors:
%       - Pedro Saa     2014 original code adapted from Qi et al. 2009

A = [];
if(isempty(B))
    A = C;
else
    NumLeft = size(B, 1);
    NumRight = size(C, 1);
    for i = 1:NumLeft
        for j = 1:NumRight
            temp = MyTimes(B(i, :), C(j, :));
            if(~isempty(temp))               
                rowNum = strmatch(round(temp),round(A),'exact');
                if(isempty(rowNum))
                    A = [A;temp];
                else
                    A(rowNum,:) = [];
                end
            end
        end
    end
end
A = sortrows(A);
    function  [T] = MyTimes(B, C)     
        if(sum(find(B==C))>0)
            T = [];
        else
            T = [B, C];
        end
        T = sort(T);
    end
end