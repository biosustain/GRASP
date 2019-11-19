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

% 1. Get output file handler
reactionName = [reactionName,num2str(strucIdx)];
currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
filepath = fullfile(currentPath{1}, '..', '..', 'temp', 'reactions', [reactionName,'.m']);

try
    fid = fopen(filepath, 'w'); 
catch
    error(['File not found: ', filepath, ...
            newline, ...
           'Please make sure the folder ', fullfile(currentPath{1}, '..', '..', 'temp'), ...
           ' exists.']) ; 
end

% 2. Write exchange mechanism
c = '%';
fprintf(fid,['function v = ',reactionName,'(SC,S,PC,P,K)\n']);
fprintf(fid,'%s Mass action definition \n',c);

rateLaw = 'v = K(1)*prod(S.^SC)-K(2)*prod(P.^PC);';

fprintf(fid,rateLaw);

fclose(fid);