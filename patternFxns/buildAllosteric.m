function buildAllosteric(metList,reactionName,negEffectors,posEffectors)
%--------------------------------------------------------------------------
% Builds allosteric reaction file
% 
% Inputs: (metList)  list with all mets concentrations showing up in the
%                    pseudo-first-order rate constants 
%    (reactionName)  reaction name
%    (negEffectors)  list with negative effectors
%    (posEffectors)  list with positive effectors
%
% Outputs:      --   writen .m file with the reaction mechanism
%--------------------- Pedro Saa 2016 -------------------------------------
% 1. Write initial parameters
try
    fid = fopen(['reactions/',reactionName,'.m'],'w'); 
catch
    fid = fopen([reactionName,'.m'],'w'); 
end
c = '%';
if (isempty(metList))
    fprintf(fid,['function v = ',reactionName,'(X,negEff,posEff,Kr,KposEff,KnegEff,L,n) \n']);
end
if(~isempty(metList))
    if ~isempty(negEffectors) && ~isempty(posEffectors)
        fprintf(fid,['function v = ',reactionName,'(X,negEff,posEff,Kr,KnegEff,KposEff,L,n) \n']);   
    elseif ~isempty(negEffectors)
        fprintf(fid,['function v = ',reactionName,'(X,negEff,Kr,KnegEff,L,n) \n']);   
    elseif ~isempty(posEffectors)
        fprintf(fid,['function v = ',reactionName,'(X,posEff,Kr,KposEff,L,n) \n']);
    else
        fprintf(fid,['function v = ',reactionName,'(X,Kr,L,n) \n']);
    end        
end
fprintf(fid,'%s Parameters definition \n',c);
strform = strcat(reactionName,'Catalytic');

% Active form reaction rate
fprintf(fid,['[vR,eR] = ',strform,'(X,Kr); \n']);

% 2. Print reaction terms
fprintf(fid,'Q = L*eR.^n; \n');
if ~isempty(negEffectors)    
    fprintf(fid,'KnegEff = KnegEff(ones(size(negEff,2),1),:); \n');
    fprintf(fid,'Q = Q.*((1 + sum(negEff''./KnegEff,2)).^n)''; \n');
end
if ~isempty(posEffectors)
    fprintf(fid,'KposEff = KposEff(ones(size(posEff,2),1),:); \n');
    fprintf(fid,'Q = Q.*((1 + sum(posEff''./KposEff,2)).^-n)''; \n');
end
fprintf(fid,'%s Reaction rate \n',c);
fprintf(fid,'v = n*vR./(1 + Q);');
fclose(fid);