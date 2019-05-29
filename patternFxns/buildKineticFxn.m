function [rxnMetLinks,freeVars,metsActive] = buildKineticFxn(ensemble,kineticFxn,strucIdx)
%--------------------------------------------------------------------------
% Builds kinetic reaction file
%
% Inputs: ensemble    ensemble (structure), kinetic fxn name, strucIdx
%
% Outputs:    -       (writen .m file with the reaction mechanism)
%-----------------Pedro Saa 2014, Marta Matos 2018-------------------------
% Define active species (mets/enzymes)
metsActive = ensemble.metsSimulated(~ismember(ensemble.metsSimulated,ensemble.metsFixed));
enzActive  = ensemble.activeRxns(~ismember(ensemble.activeRxns,ensemble.kinInactRxns));
totalEvals = numel(metsActive) + numel(enzActive) + 1;
freeVars   = [ensemble.mets(metsActive);ensemble.rxns(enzActive)];                             % return indexes of the free variables

% Write initial parameters
c = '%';
fid = fopen(['reactions',num2str(strucIdx),'/',kineticFxn,'.m'],'w');
fprintf(fid,['function [f,grad] = ',kineticFxn,'(x,model,fixedExch,Sred,kinInactRxns,subunits,flag)\n']);
fprintf(fid,'%s Pre-allocation of memory\n',c);
fprintf(fid,'h = 1e-8;\n');									% Step length
fprintf(fid,'%s Defining metabolite and enzyme species\n',c);
fprintf(fid,'if flag==1\n');
fprintf(fid,'x = x(:);\n');									% Column vector
fprintf(fid,['v = zeros(',num2str(size(ensemble.Sred,2)),',',num2str(totalEvals),');\n']);      % Preallocation of memory (rxns)
fprintf(fid,['E = zeros(',num2str(size(ensemble.Sred,2)),',',num2str(totalEvals),');\n']);      % Preallocation of memory (enz)
fprintf(fid,['x = [x,x(:,ones(1,',num2str(totalEvals-1),')) + diag(h*1i*ones(',num2str(totalEvals-1),',1))];\n']);      % Preallocation of memory (free vars)
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
if ~isempty(ensemble.kinInactRxns)
    fprintf(fid,['E(kinInactRxns,:) = fixedExch(:,ones(1,size(x,2)));\n']);  % Define fixed protein concentrations
end

% Define equations and determine connectivity between reactions and
% metabolites
rxnMetLinks{length(ensemble.activeRxns)} = [];
k = 1;
fprintf(fid,'%s Reaction rates\n',c);
for i = 1:numel(ensemble.activeRxns)

    % If binding/release order is provided
    % Does not take into account different models
    if ~isempty(ensemble.subOrder{strucIdx}{i})||~isempty(ensemble.prodOrder{strucIdx}{i})
        reactants = [];
        w = sum(length(ensemble.inhibitors{1,1}{i})+length(ensemble.activators{1,1}{i})+length(ensemble.subOrder{1,1}{i})+length(ensemble.prodOrder{1,1}{i}));
        elemCount = 1;
        prodCount = 1;
        subCount = 1;
        inhCount = 1;
        actCount = 1;
        for react = 1:length(ensemble.metLists{i,1})
            if ~isempty(strfind(ensemble.metLists{i,1}{react}, 'A'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'B'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'C'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'D'))
                if elemCount < w
                    reactants = [reactants, ensemble.subOrder{1, 1}{i}{subCount}, ';'];
                    subCount = subCount+1;
                    elemCount = elemCount+1;
                else
                    reactants = [reactants, ensemble.subOrder{1, 1}{i}{subCount}];
                end
            elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'P'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'Q'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'R'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'S'))
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

        for reactantI = 1:size(reactants, 2)
            if ismember(reactants{reactantI},ensemble.mets(ensemble.metsFixed))
                reactants{reactantI} = strcat('ones(1,size(x,2))');
            end
        end


        if strcmp('massAction',ensemble.rxnMechanisms{strucIdx}(i))
            reactants = [strjoin(reactants, ',') ','];
        else
            reactants = strjoin(reactants, ';');
        end
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
        substrates  = substrates(ismember(substrates,ensemble.mets(ensemble.metsSimulated)));         % Extract only active substrates
        stoicCoeffs = abs(ensemble.S(ismember(ensemble.mets,substrates),ensemble.activeRxns(i)));

        % Substrates: Check the stoic coeff (relevant only if greater than 1)
        if any(stoicCoeffs>1)
            for w = 1:numel(substrates)
                if stoicCoeffs(w)>1
                    for u = 2:stoicCoeffs(w)
                        substrates = [substrates;substrates(w)];
                    end
                end
            end
        end

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
        products    = products(ismember(products,ensemble.mets(ensemble.metsSimulated)));               % Extract only active products
        stoicCoeffs = abs(ensemble.S(ismember(ensemble.mets,products),ensemble.activeRxns(i)));


        % Products: Check the stoic coeff (relevant only if greater than 1)
        if any(stoicCoeffs>1)
            for w = 1:numel(products)
                if stoicCoeffs(w)>1
                    for u = 2:stoicCoeffs(w)
                        products = [products;products(w)];
                    end
                end
            end
        end

        % Non-enzymatic reactions (diffusion)
        if strcmp('diffusion',ensemble.rxnMechanisms{strucIdx}(i))
            if ~isempty(substrates)
                reactants = substrates{1};
            else
                reactants = products{1};
            end
            rxnMetLinks{i} = reactants;                                                           % There is a link if rxn is diffusion

        elseif strcmp('massAction',ensemble.rxnMechanisms{strucIdx}(i))
            if ~ismember(substrates{1},ensemble.mets(ensemble.metsFixed))
                reactants      = [reactants,substrates{1},','];
                rxnMetLinks{i} = [rxnMetLinks{i},substrates(1)];
            else
                reactants = [reactants,'ones(1,size(x,2)),'];
            end
            if ~ismember(products{1},ensemble.mets(ensemble.metsFixed))
                reactants      = [reactants,products{1},','];
                rxnMetLinks{i} = [rxnMetLinks{i},products(1)];
            else
                reactants = [reactants,'ones(1,size(x,2)),'];
            end

        % Enzymatic reactions
        elseif  ~strcmp('fixedExchange',ensemble.rxnMechanisms{strucIdx}(i))

            w = sum(length(ensemble.inhibitors{1,1}{i})+length(ensemble.activators{1,1}{i})+length(substrates)+length(products));
            elemCount = 1;
            prodCount = 1;
            subCount = 1;
            inhCount = 1;
            actCount = 1;
            for react = 1:length(ensemble.metLists{i,1})
                if ~isempty(strfind(ensemble.metLists{i,1}{react}, 'A'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'B'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'C'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'D'))
                    if elemCount < w
                        reactants = [reactants,substrates{subCount}, ';'];
                        rxnMetLinks{i} = [rxnMetLinks{i}, substrates{subCount}];
                        subCount = subCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants, substrates{subCount}];
                        rxnMetLinks{i} = [rxnMetLinks{i}, substrates{subCount}];
                    end
                elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'P'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'Q'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'R'))||~isempty(strfind(ensemble.metLists{i,1}{react}, 'S'))
                    if elemCount < w
                        reactants = [reactants,products{prodCount},';'];
                        rxnMetLinks{i} = [rxnMetLinks{i}, products{prodCount}];
                        prodCount = prodCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants,products{prodCount}];
                        rxnMetLinks{i} = [rxnMetLinks{i}, products{prodCount}];
                    end
                elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'I'))
                    if elemCount < w
                        reactants = [reactants, ensemble.inhibitors{1, 1}{i}{inhCount},';'];
                        rxnMetLinks{i} = [rxnMetLinks{i}, ensemble.inhibitors{1, 1}{i}{inhCount}];
                        inhCount = inhCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants, ensemble.inhibitors{1, 1}{i}{inhCount}];
                        rxnMetLinks{i} = [rxnMetLinks{i}, ensemble.inhibitors{1, 1}{i}{inhCount}];
                    end
                elseif ~isempty(strfind(ensemble.metLists{i,1}{react}, 'Z'))
                    if elemCount < w
                        reactants = [reactants,ensemble.activators{1, 1}{i}{actCount}, ';'];
                        rxnMetLinks{i} = [rxnMetLinks{i}, ensemble.activators{1, 1}{i}{actCount}];
                        actCount = actCount+1;
                        elemCount = elemCount+1;
                    else
                        reactants = [reactants,ensemble.activators{1, 1}{i}{actCount}];
                        rxnMetLinks{i} = [rxnMetLinks{i}, ensemble.activators{1, 1}{i}{actCount}];
                    end
                end
            end

        reactants = strsplit(reactants, ';');

        for reactantI = 1:size(reactants, 2)
            if ismember(reactants{reactantI},ensemble.mets(ensemble.metsFixed))
                reactants{reactantI} = strcat('ones(1,size(x,2))');
            end
        end

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
            fprintf(fid,['v(%i,:) = ',ensemble.rxns{i},num2str(strucIdx),'([],[]);\n'],i);
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

    % Add enzyme links to the rxnMet links (only if the reaction is not an exch)
    if ~isempty(reactants)
        rxnMetLinks{i} = [rxnMetLinks{i},ensemble.rxns{ensemble.activeRxns(i)}];
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
