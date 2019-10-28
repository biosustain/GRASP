function buildMassAction(reactionName,strucIdx)
% Build mass action kinetic function.
% Note that this rate law assumes only one substrate and one product.
%
%
% USAGE:
%
%    buildMassAction(reactionName, strucIdx)
%
% INPUT:
%    reactionName (char):   reaction name
%    strucIdx (int):        ID of the model structure
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
fprintf(fid,['function v = ',reactionName,'(S,P,K)\n']);
fprintf(fid,'%s Mass action definition \n',c);
fprintf(fid,'v = K(1)*prod(S,1)-K(2)*prod(P,1);');
fclose(fid);