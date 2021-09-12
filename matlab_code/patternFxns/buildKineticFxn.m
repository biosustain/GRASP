function freeVars = buildKineticFxn(ensemble,kineticFxn,strucIdx)
% Builds kinetic model .m file that is used to actually run the model.
%
%
% USAGE:
%
%    freeVars = buildKineticFxn(ensemble, kineticFxn, strucIdx)
%
% INPUT:
%    ensemble (struct):   model ensemble, see buildEnsemble for fields description
%    kineticFxn (char):	  name of the kinetic function
%    strucIdx (int):      ID of the model structure
%
% OUTPUT:
%    freeVars (char cell):      [TODO Pedro]
%    written .m file with the model
%
% .. Authors:
%       - Pedro Saa         2016 original code 
%       - Marta Matos       2018 extended for promiscuous reactions
%       - Nicholas Cowie	2019 extended for isoenzymes

% Define active species (mets/enzymes)
metsActive = ensemble.metsActive;
enzActive  = ensemble.activeRxns;
%enzActive  = ensemble.activeRxns(~ismember(ensemble.activeRxns,ensemble.kinInactRxns));
totalEvals = numel(metsActive) + numel(enzActive) + 1;
freeVars   = [ensemble.mets(metsActive);ensemble.rxns(enzActive)];                             % return indexes of the free variables


% Get output file handler
currentPath = regexp(mfilename('fullpath'), '(.*)[/\\\\]', 'match');
filepath = fullfile(currentPath{1}, '..', '..', 'reactions', [ensemble.description, '_', num2str(strucIdx)], [kineticFxn,'.m']);

fid = fopen(filepath, 'w'); 
if fid == -1
    error(['File not found: ', filepath, ...
           newline, ...
           'Please make sure the folder ', fullfile(currentPath{1}, '..', '..', 'reactions'), ...
           ' exists.']) 
end

% Write initial parameters
c = '%';
fprintf(fid,['function [f,grad] = ',kineticFxn,'(x,xconst,model,fixedExch,Sred,kinInactRxns,subunits,flag)\n']);
fprintf(fid,'%s Pre-allocation of memory\n',c);
fprintf(fid,'h = 1e-8;\n');									% Step length
fprintf(fid,'%s Defining metabolite and enzyme species\n',c);
fprintf(fid,'if flag==1\n');
fprintf(fid,'x = x(:);\n');									% Column vector
fprintf(fid,'xconst = xconst(:);\n');									% Column vector
fprintf(fid,['v = zeros(',num2str(size(ensemble.Sred,2)),',',num2str(totalEvals),');\n']);      % Preallocation of memory (rxns)
fprintf(fid,['E = zeros(',num2str(size(ensemble.Sred,2)),',',num2str(totalEvals),');\n']);      % Preallocation of memory (enz)
fprintf(fid,['x = [x,x(:,ones(1,',num2str(totalEvals-1),')) + diag(h*1i*ones(',num2str(totalEvals-1),',1))];\n']);      % Preallocation of memory (free vars)
fprintf(fid,['xconst = [xconst,xconst(:,ones(1,',num2str(totalEvals-1),'))];\n']);      % Preallocation of memory (constant metabolites)
fprintf(fid,'else\n');
fprintf(fid,['v = zeros(',num2str(size(ensemble.Sred,2)),',size(x,2));\n']);      % Preallocation of memory (rxns)
fprintf(fid,['E = zeros(',num2str(size(ensemble.Sred,2)),',size(x,2));\n']);      % Preallocation of memory (enz)
fprintf(fid,'end\n');

% Define metabolite species
fprintf(fid,'%s Defining metabolite and enzyme species\n',c);

i = 1; r = 1;
for j = 1:numel(metsActive)+numel(enzActive)
    if j<=numel(metsActive)
        fprintf(fid,[ensemble.mets{metsActive(j)},' = x(%i,:);\n'],i);                     % Only simulated metabolites in the kinetic fxn
    else
        fprintf(fid,['E(',num2str(enzActive(r)),',:) = x(%i,:);\n'],i);                    % Only simulated enzymes in the kinetic fxn
        r = r+1;
    end
    i = i+1;
end

for i = 1:numel(ensemble.metsFixed)
    fprintf(fid,[ensemble.mets{ensemble.metsFixed(i)},' = xconst(%i,:);\n'],i);                     % Only simulated metabolites in the kinetic fxn
end


fixedExchangeI  = 1;  % To keep track of the fixed exchange index

% Define equations and determine connectivity between reactions and
% metabolites
k = 1;
fprintf(fid,'%s Reaction rates\n',c);
for i = 1:numel(ensemble.activeRxns)

    % If binding/release order is provided
    % Does not take into account different models
    if (~isempty(ensemble.subOrder{strucIdx}{i}) || ~isempty(ensemble.prodOrder{strucIdx}{i})) && sum(~ismember({'massAction', 'fixedExchange', 'freeExchange', 'diffusion'}, ensemble.rxnMechanisms{strucIdx}{i})) == 4 
        reactants = [];
        w = sum(length(ensemble.inhibitors{1,1}{i})+length(ensemble.activators{1,1}{i})+length(ensemble.subOrder{1,1}{i})+length(ensemble.prodOrder{1,1}{i}));
        elemCount = 1;
        prodCount = 1;
        subCount = 1;
        inhCount = 1;
        actCount = 1;
        for react = 1:length(ensemble.metLists{i,1})
            if (~isempty(strfind(ensemble.metLists{i,1}{react}, 'A')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'B')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'C')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'D')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'E')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'F')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'G')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'H')) || ...
                ~isempty(strfind(ensemble.metLists{i,1}{react}, 'J')))
            
                if elemCount < w
                    reactants = [reactants, ensemble.subOrder{1, 1}{i}{subCount}, ';'];
                    subCount = subCount+1;
                    elemCount = elemCount+1;
                else
                    reactants = [reactants, ensemble.subOrder{1, 1}{i}{subCount}];
                end
            elseif (~isempty(strfind(ensemble.metLists{i,1}{react}, 'P')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'Q')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'R')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'S')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'T')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'U')) || ...    
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'V')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'W')) || ...
                    ~isempty(strfind(ensemble.metLists{i,1}{react}, 'X')))
                
                if elemCount < w
                    reactants = [reactants,ensemble.prodOrder{1, 1}{i}{prodCount},';'];
                    prodCount = prodCount+1;
                    elemCount = elemCount+1;
                else
                    reactants = [reactants,ensemble.prodOrder{1, 1}{i}{prodCount}];
                end
            elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'I'))
                if elemCount < w
                    reactants = [reactants, ensemble.inhibitors{1, 1}{i}{inhCount},';'];
                    inhCount = inhCount+1;
                    elemCount = elemCount+1;
                else
                    reactants = [reactants, ensemble.inhibitors{1, 1}{i}{inhCount}];
                end
            elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'Z'))
                if elemCount < w
                    reactants = [reactants,ensemble.activators{1, 1}{i}{actCount}, ';'];
                    actCount = actCount+1;
                    elemCount = elemCount+1;
                else
                    reactants = [reactants,ensemble.activators{1, 1}{i}{actCount}];
                end
            end
        end
        
        reactants = strsplit(reactants, ';');
        reactants = strjoin(reactants, ';');
        
    else
        
        reactants  = [];

        % Extract and organize substrates
        metsInd = 1:size(ensemble.mets);
        if size(ensemble.promiscuity{strucIdx}{i}) > 0
            substratesInd = [];
            for rxnI = ensemble.promiscuity{strucIdx}{i}
                substratesIndTemp = metsInd(ensemble.S(:,ensemble.activeRxns(rxnI))<0);
                substratesInd = [substratesInd substratesIndTemp];
            end
            substrates = ensemble.mets(substratesInd);
        else
            substrates  = (ensemble.mets(ensemble.S(:,ensemble.activeRxns(i))<0));
        end
        
        stoicCoeffsSub = abs(ensemble.S(ismember(ensemble.mets,substrates),ensemble.activeRxns(i)));
      
        % Extract and organize products
        if size(ensemble.promiscuity{strucIdx}{i}) > 0
            productsInd = [];
            for rxnI = ensemble.promiscuity{strucIdx}{i}
                productsIndTemp = (metsInd(ensemble.S(:,ensemble.activeRxns(rxnI))>0));
                productsInd = [productsInd productsIndTemp];
            end
            products = ensemble.mets(productsInd);
        else
            products    = (ensemble.mets(ensemble.S(:,ensemble.activeRxns(i))>0));
        end
        
        stoicCoeffsProd = abs(ensemble.S(ismember(ensemble.mets,products),ensemble.activeRxns(i)));

        % Non-enzymatic reactions (diffusion)
        if strcmp('diffusion',ensemble.rxnMechanisms{strucIdx}(i))
            if ~isempty(substrates)
                reactants = substrates{1};
            else
                reactants = products{1};
            end

        elseif strcmp('massAction',ensemble.rxnMechanisms{strucIdx}(i))
                    
            subList = [];
            subCoefList = [];
            for subI=1:numel(substrates)
                subCoefList = [subCoefList, num2str(stoicCoeffsSub(subI)), '*ones(1,size(x,2));'];
                subList     = [subList, substrates{subI}, ';'];
            end
            
            prodList = [];
            prodCoefList = [];
            for prodI=1:numel(products)
                prodCoefList = [prodCoefList, num2str(stoicCoeffsProd(prodI)), '*ones(1,size(x,2));'];
                prodList     = [prodList, products{prodI}, ';'];
            end
            
            
            reactants = ['[', subCoefList(1:end-1), ...
                         '],[', subList(1:end-1), '],[',  ...
                         prodCoefList(1:end-1), ...
                         '],[',prodList(1:end-1), '],'];
            
            
        % Enzymatic reactions
        elseif  ~strcmp('fixedExchange',ensemble.rxnMechanisms{strucIdx}(i))

            w = sum(length(ensemble.inhibitors{1,1}{i})+length(ensemble.activators{1,1}{i})+length(substrates)+length(products));
            elemCount = 1;
            prodCount = 1;
            subCount = 1;
            inhCount = 1;
            actCount = 1;
            for react = 1:length(ensemble.metLists{i,1})
                 if (~isempty(strfind(ensemble.metLists{i,1}{react}, 'A')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'B')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'C')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'D')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'E')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'F')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'G')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'H')) || ...
                     ~isempty(strfind(ensemble.metLists{i,1}{react}, 'J')))
                
                    if elemCount < w
                        reactants = [reactants,substrates{subCount}, ';'];
                        subCount = subCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants, substrates{subCount}];
                    end
                elseif (~isempty(strfind(ensemble.metLists{i,1}{react}, 'P')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'Q')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'R')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'S')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'T')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'U')) || ...    
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'V')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'W')) || ...
                        ~isempty(strfind(ensemble.metLists{i,1}{react}, 'X')))
                
                    if elemCount < w
                        reactants = [reactants,products{prodCount},';'];
                        prodCount = prodCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants,products{prodCount}];
                    end
                elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'I'))
                    if elemCount < w
                        reactants = [reactants, ensemble.inhibitors{1, 1}{i}{inhCount},';'];
                        inhCount = inhCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants, ensemble.inhibitors{1, 1}{i}{inhCount}];
                    end
                elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'Z'))
                    if elemCount < w
                        reactants = [reactants,ensemble.activators{1, 1}{i}{actCount}, ';'];
                        actCount = actCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants,ensemble.activators{1, 1}{i}{actCount}];
                    end
                end
            end

        reactants = strsplit(reactants, ';');
        reactants = strjoin(reactants, ';');
        end
    end

    % Allosteric reaction
    if ensemble.allosteric{strucIdx}(i)
        negEffectors = [];
        posEffectors = [];

        % Both positive and negative effectors
        if ~isempty(ensemble.negEffectors{strucIdx}{i}) && ~isempty(ensemble.posEffectors{strucIdx}{i})
            for j = 1:length(ensemble.negEffectors{strucIdx}{i})-1
                negEffectors = [negEffectors,char(ensemble.negEffectors{strucIdx}{i}(j)),';'];
            end
            negEffectors = [negEffectors,char(ensemble.negEffectors{strucIdx}{i}(length(ensemble.negEffectors{strucIdx}{i})))];
            for j = 1:length(ensemble.posEffectors{strucIdx}{i})-1
                posEffectors = [posEffectors,char(ensemble.posEffectors{strucIdx}{i}(j)),';'];
            end
            posEffectors = [posEffectors,char(ensemble.posEffectors{strucIdx}{i}(length(ensemble.posEffectors{strucIdx}{i})))];
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'([',reactants,'],[',negEffectors,'],[',posEffectors,'],model.rxnParams(%i).kineticParams,model.rxnParams(%i).KnegEff,model.rxnParams(%i).KposEff,model.rxnParams(%i).L,subunits(%i));\n'],i,k,k,k,k,k);
            k = k+1;

        % Only negative effectors
        elseif ~isempty(ensemble.negEffectors{strucIdx}{i})
            for j = 1:length(ensemble.negEffectors{strucIdx}{i})-1
                negEffectors = [negEffectors,char(ensemble.negEffectors{strucIdx}{i}(j)),';'];
            end
            negEffectors = [negEffectors,char(ensemble.negEffectors{strucIdx}{i}(length(ensemble.negEffectors{strucIdx}{i})))];
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'([',reactants,'],[',negEffectors,'],model.rxnParams(%i).kineticParams,model.rxnParams(%i).KnegEff,model.rxnParams(%i).L,subunits(%i));\n'],i,k,k,k,k);
            k = k+1;

        % Only positive effectors
        elseif ~isempty(ensemble.posEffectors{strucIdx}{i})
            for j = 1:length(ensemble.posEffectors{strucIdx}{i})-1
                posEffectors = [posEffectors,char(ensemble.posEffectors{strucIdx}{i}(j)),';'];
            end
            posEffectors = [posEffectors,char(ensemble.posEffectors{strucIdx}{i}(length(ensemble.posEffectors{strucIdx}{i})))];
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'([',reactants,'],[',posEffectors,'],model.rxnParams(%i).kineticParams,model.rxnParams(%i).KposEff,model.rxnParams(%i).L,subunits(%i));\n'],i,k,k,k,k);
            k = k+1;

        % No effectors, only cooperative binding
        else
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'([',reactants,'],model.rxnParams(%i).kineticParams,model.rxnParams(%i).L,subunits(%i));\n'],i,k,k,k);
            k = k+1;
        end

        % Non-allosteric reaction
    else
        if strcmp('fixedExchange',ensemble.rxnMechanisms{strucIdx}(i))
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx), '(fixedExch(%i), size(x,2));\n'], i, fixedExchangeI);
            fixedExchangeI = fixedExchangeI + 1;
        elseif strcmp('freeExchange',ensemble.rxnMechanisms{strucIdx}(i))
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'(model.rxnParams(%i).kineticParams,',num2str(totalEvals),');\n'],i,k);
            k = k+1;
		elseif strcmp('massAction',ensemble.rxnMechanisms{strucIdx}(i))
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'(',reactants,'model.rxnParams(%i).kineticParams);\n'],i,k);
            k = k+1;
        else
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'([',reactants,'],model.rxnParams(%i).kineticParams);\n'],i,k);
            k = k+1;
        end
    end

end

% Definition of final rates
fprintf(fid,'if flag==1\n');
fprintf(fid,'%s Final rates\n',c);
fprintf(fid,'y = sum((Sred*(E.*v)).^2);\n');
fprintf(fid,'f = real(y(1));\n');
fprintf(fid,'if (nargout>1) %s gradient is required\n',c);
fprintf(fid,'grad = imag(y(2:end))/h;\n');
fprintf(fid,'end\n');
fprintf(fid,'else\n');
fprintf(fid,'f = E.*v;\n');
fprintf(fid,'grad = [];\n');
fprintf(fid,'end');
fclose(fid);
