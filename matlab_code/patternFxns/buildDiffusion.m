function buildDiffusion(reactionName,strucIdx)
%--------------------------------------------------------------------------
% Build diffusion kinetic fxn
%
% Inputs: ensemble    (ensemble structure)
%
% Outputs:    -       (writen .m file with the reaction mechanism)
%------------------------Pedro Saa 2016------------------------------------
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
fprintf(fid,['function v = ',reactionName,'(X,K)\n']);
fprintf(fid,'%s Difusion reaction definition \n',c);
fprintf(fid,'v = K*X;\n');
fclose(fid);