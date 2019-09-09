function buildFreeExchange(filename,strucIdx)
%--------------------------------------------------------------------------
% Build free exchange kinetic fxn
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
fprintf(fid,['function v = ',filename,'(K,numConditions)\n']);
fprintf(fid,'%s Free exchange definition \n',c);
fprintf(fid,'v = K(ones(1,numConditions));\n');
fclose(fid);