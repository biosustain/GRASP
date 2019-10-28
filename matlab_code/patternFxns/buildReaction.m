function buildReaction(state,rateList,metList,numTerm,prodNum,reactionName, promiscuousRxnI)
% Builds rate reaction file.
%
% Based on https://doi.org/10.1186/1471-2105-10-238 
%
%
% USAGE:
%
%    buildReaction(state, rateList, metList, numTerm, prodNum, reactionName, promiscuousRxnI)
%
% INPUT:
%    state (char cell):       output from reactionPatterm     
%    rateList (char cell):	  list with all the rate constants in the system
%    metList (char cell):     list with all mets concentrations showing up in thepseudo-first-order rate constants
%    numTerm (char cell):     numerator terms
%    prodNum (int vector):	  products terms
%    reactionName (char):     reaction name
%    promiscuousRxnI (int):	  index of promiscuous reaction
%
% OUTPUT:
%    written .m file with the reaction mechanism
%
% .. Authors:
%       - Pedro Saa     2016 original code  adapted from Qi et al
%       - Marta Matos   2018 extended for promiscuous reactions

% 1. Write initial parameters
try
    currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
    filepath = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions', [reactionName,'.m']);
    fid = fopen(filepath, 'w'); 
catch
    fid = fopen([reactionName,'.m'],'w'); 
end
c = '%';
if (isempty(metList))
    fprintf(fid,['function [v,E1,E2] = ',reactionName,'(K)\n']);
end
if(~isempty(metList))
    fprintf(fid,['function [v,E1,E2] = ',reactionName,'(X,K)\n']);   
    fprintf(fid,'%s Metabolites definition \n',c);
    len = length(metList);
    for i = 1:len
        temp = metList{i};
        temp = regexp(temp,'*','split');
        fprintf(fid,'%s = X(%i,:);\n',temp{2},i);
    end
end
fprintf(fid,'%s Parameters definition K \n', c);
len = length(rateList);
for i = 1:len
    temp = rateList{i};
    fprintf(fid,'%s = K(%i);\n',temp,i);
end
fprintf(fid,'%s  Numerator terms\n',c);
len = length(state);
for i = 1:len
    temp = char(state{i});
    loc = strfind(temp,'=');
    fprintf(fid,'E%i =%s\n',i,temp(loc+1:end));
end

% 2. Print rate terms
fprintf(fid,'%s Denominator terms\n',c);
fprintf(fid,'D = E1');
for i = 2:len-1
    fprintf(fid,'+E%i',i);
end
fprintf(fid,'+E%s;\n',num2str(len));
fprintf(fid,'%s Enzyme abundances terms\n',c);
for i = 1:len
    fprintf(fid,'E%i = E%i./D;\n',i,i);
end
fprintf(fid,'%s Reaction rate \n', c);
fprintf(fid,'v = ');

for i = 1:size(prodNum,1)
    if promiscuousRxnI == 0 || promiscuousRxnI == i
        fprintf(fid,numTerm{i,1});
        fprintf(fid,'.*E%i',prodNum(i,1));
        fprintf(fid,numTerm{i,2});
        fprintf(fid,'.*E%i',prodNum(i,2));
    end
end
fprintf(fid,';');
fclose(fid);