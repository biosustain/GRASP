function buildReaction(state,rateList,metList,numTerm,prodNum,filename, promiscuousRxnI)
%--------------------------------------------------------------------------
% Builds rate reaction file
% 
% Inputs:   (state)  output from reactionPatterm         
%        (rateList)  list with all the rate constants in the system
%         (metList)  list with all mets concentrations showing up  in the
%                    pseudo-first-order rate constants 
%         (numTerm)  numerator terms
%         (prodNum)  products terms
%        (filename)  pattern name
%
% Outputs:      --   writen .m file with the reaction mechanism
%--------Pedro Saa 2016, adapted from Qi et al. 2009, Marta Matos 2018-----
% 1. Write initial parameters
try
    fid = fopen(['reactions/',filename,'.m'],'w'); 
catch
    fid = fopen([filename,'.m'],'w'); 
end
c = '%';
if (isempty(metList))
    fprintf(fid,['function [v,E1,E2] = ',filename,'(K)\n']);
end
if(~isempty(metList))
    fprintf(fid,['function [v,E1,E2] = ',filename,'(X,K)\n']);   
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