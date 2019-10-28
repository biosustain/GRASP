function cellWithVariables = fixVariableNames(cellWithVariables, prefix, flag)
% Takes in a cell where the first column are metabolite or reaction names
% and adds the specified prefix to each name.
%
% For instance, if a metabolite is named nad and the prefix is specified 
% as 'm', then the metabolite new name will be m_nad" ?
%
% If *flag* is defined (independently of its value) the same is done on the 
% kinetics sheet for the columns order, promiscuous, inhibitors, 
% activators, negative effector, positive effector, allosteric.
%
%
% USAGE:
%
%    cellWithVariables = fixVariableNames(cellWithVariables, prefix, flag)
%
% INPUT:
%    cellWithVariables (char cell):  variable names
%    prefix (char):                  prefix to add to variable names
%    flag (logical):                 whether or not to do this for the kinetics sheet.
%                                  
%
% OUTPUT:
%    cellWithVariables (cell):	variable names with prefixes
%
% .. Authors:
%       - Marta Matos	2019 original code


for row=2:size(cellWithVariables,1)
    cellWithVariables{row,1} = strcat(prefix, '_', cellWithVariables{row,1});
end

% add m and r to other columns in the kinetics sheet
if (nargin > 2)
    for col=3:(size(cellWithVariables,2)-3)

        if cellWithVariables{1,col} == "promiscuous"
            prefix = "r";
        else
            prefix = "m";
        end
        
        for row=2:size(cellWithVariables,1)            
            if (cellWithVariables{row,col} ~= "")
                
            	char_list = strsplit(cellWithVariables{row,col}, ' ');
                new_entry = strcat(prefix, '_', char_list{1}); 
                
                for i=2:size(char_list,2)
                    new_entry = strcat(new_entry, {' '}, prefix, '_', char_list{i});
                end
                cellWithVariables{row,col} = char(new_entry);
            end
        end     
    end
end

end

