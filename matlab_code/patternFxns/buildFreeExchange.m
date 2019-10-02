function buildFreeExchange(reactionName,strucIdx)
% Build free exchange kinetic fxn
%
% USAGE:
%
%    buildFreeExchange(reactionName, strucIdx)
%
% INPUTS:
%    reactionName (`char`):   reaction name
%    strucIdx (`int`):        ID of the model structure
%
% OUTPUT:
%    written .m file with the reaction mechanism
%
% .. Authors:
%       - Pedro Saa     2016 original code 

reactionName = [reactionName,num2str(strucIdx)];
try
    currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
    filepath = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions', [reactionName,'.m']);
    fid = fopen(filepath, 'w'); 
catch
    fid = fopen([reactionName,'.m'],'w'); 
end
%% 1. Write exchange mechanism
c = '%';
fprintf(fid,['function v = ',reactionName,'(K,numConditions)\n']);
fprintf(fid,'%s Free exchange definition \n',c);
fprintf(fid,'v = K(ones(1,numConditions));\n');
fclose(fid);