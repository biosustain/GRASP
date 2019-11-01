function buildDiffusion(reactionName, strucIdx)
% Build diffusion kinetic function.
%
%
% USAGE:
%
%    buildDiffusion(reactionName, strucIdx)
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

% 1. Get output file handler
reactionName = [reactionName,num2str(strucIdx)];
currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
filepath = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions', [reactionName,'.m']);

try
    fid = fopen(filepath, 'w'); 
catch
    error(['File not found: ', filepath, ...
            newline, ...
           'Please make sure the folder ', fullfile(currentPath{1}, '..', '..', 'temp'), ...
           ' exists.'])
end

% 2. Write exchange mechanism
c = '%';
fprintf(fid,['function v = ',reactionName,'(X,K)\n']);
fprintf(fid,'%s Difusion reaction definition \n',c);
fprintf(fid,'v = K*X;\n');
fclose(fid);