function buildMassAction(filename,strucIdx)
%--------------------------------------------------------------------------
% Build mass action kinetic fxn
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
fprintf(fid,['function v = ',filename,'(S,P,K)\n']);
fprintf(fid,'%s Mass action definition \n',c);
fprintf(fid,'v = K(1)*prod(S,1)-K(2)*prod(P,1);');
fclose(fid);