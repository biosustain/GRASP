function buildMassAction(reactionName,substrateList,productList,strucIdx)
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
fprintf(fid,['function v = ',reactionName,'(S,P,K)\n']);
fprintf(fid,'%s Mass action definition \n',c);

if ~isempty(substrateList)
    rateLaw = 'v = K(1)';

    for subI=1:size(substrateList,2)
        if contains(substrateList{subI}, '*')
            coefMet = strsplit(substrateList{subI}, '*');
            rateLaw = strcat(rateLaw, '*S(', num2str(subI), ',:)^', coefMet{1});
        else
            rateLaw = strcat(rateLaw, '*S(', num2str(subI), ',:)');
        end
    end

    rateLaw = strcat(rateLaw, '-K(2)');

    for prodI=1:size(productList,2)
        if contains(productList{prodI}, '*')
            coefMet = strsplit(productList{prodI}, '*');
            rateLaw = strcat(rateLaw, '*P(', num2str(prodI), ',:)^', coefMet{1});
        else
            rateLaw = strcat(rateLaw, '*P(', num2str(prodI), ',:)');
        end
    end

    rateLaw = strcat(rateLaw, ';');
    
else
    rateLaw = 'v = K(1)*prod(S,1)-K(2)*prod(P,1);';
end

fprintf(fid,rateLaw);

fclose(fid);