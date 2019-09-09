function buildDifusion(filename,strucIdx)
%--------------------------------------------------------------------------
% Build diffusion kinetic fxn
%
% Inputs: ensemble    (ensemble structure)
%
% Outputs:    -       (writen .m file with the reaction mechanism)
%------------------------Pedro Saa 2016------------------------------------
filename = [filename,num2str(strucIdx)];
try
    fid = fopen(['reactions/',filename,'.m'],'w'); 
catch
    fid = fopen([filename,'.m'],'w'); 
end
%% 1. Write exchange mechanism
c = '%';
fprintf(fid,['function v = ',filename,'(X,K)\n']);
fprintf(fid,'%s Difusion reaction definition \n',c);
fprintf(fid,'v = K*X;\n');
fclose(fid);