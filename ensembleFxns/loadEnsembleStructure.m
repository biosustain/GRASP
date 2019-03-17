function ensemble = loadEnsembleStructure(xlsxFile)
%--------------------------------------------------------------------------
% This function builds the data structure of the ensemble
%
% Inputs:       excel file with the invariant information of the ensemble
%
% Outputs:      ensemble basic data structure
%--------------------- Pedro Saa 2016 -------------------------------------
if nargin<1; disp('Not enough input arguments'); end                    % Check inputs

%% 1. Load general data contained in .xlsx file
xlsxFile            = [xlsxFile,'.xlsx'];                               % add extension to the file
[xData,strData]     = xlsread(xlsxFile,'general');                      % load general info
Sfull               = xlsread(xlsxFile,'stoic');                        % load full stoichiometry
[xRxns,rxnsList]    = xlsread(xlsxFile,'rxns');                         % load rxn info
[xMets,metsList]    = xlsread(xlsxFile,'mets');                         % load mets info
splitRatios         = xlsread(xlsxFile,'splitRatios');                  % load split ratios
poolConstraints     = xlsread(xlsxFile,'poolConst');                    % load pool constraints
[measRates,idxMeas] = xlsread(xlsxFile,'measRates');                    % load meas rates data
xDG_std             = xlsread(xlsxFile,'thermoRxns');                   % load thermodynamic data (rxns)
xMetsThermo         = xlsread(xlsxFile,'thermoMets');                   % load thermodynamic data (mets)
[protData,idxProt]  = xlsread(xlsxFile,'protData');                     % load expression data
[metsData,idxMets]  = xlsread(xlsxFile,'metsData');                     % load metabolite data (mets)
ineqConstraints     = xlsread(xlsxFile,'thermo_ineq_constraints');      % load ineq. thermodynamic constraints

% Add m before any metabolite name and r before any reaction name to avoid
% variables starting with a number
rxnsList = fixVariableNames(rxnsList, 'r');
metsList = fixVariableNames(metsList, 'm');
idxMeas = fixVariableNames(idxMeas, 'r');
idxProt = fixVariableNames(idxProt, 'r');
idxMets = fixVariableNames(idxMets, 'm');


% Build initial ensemble structure
ensemble.description   = strData{2,2};
ensemble.sampler       = strData{3,2};
ensemble.solver        = strData{4,2};
ensemble.numConditions = xData(1);
ensemble.numStruct     = xData(2);
ensemble.numParticles  = xData(3);
ensemble.parallel      = xData(4);
ensemble.numCores      = xData(5);
ensemble.alphaAlive    = xData(6);
robustFluxes           = xData(7);
computeThermo          = xData(8);
if isnan(xData(9))
    ensemble.tolerance = [Inf,xData(10)];
else
    ensemble.tolerance = [xData(9),xData(10)];
end
ensemble.S             = Sfull';
ensemble.rxns          = rxnsList(2:end,1);
ensemble.rxnNames      = rxnsList(2:end,2);
ensemble.exchRxns      = find(xRxns(:,1));
ensemble.activeRxns    = find(xRxns(:,2));
ensemble.isoenzymes    = rxnsList([2:end],5);
ensemble.uniqueIso       = unique(ensemble.isoenzymes(~cellfun(@isempty, ensemble.isoenzymes)));
ensemble.mets          = metsList(2:end,1);
ensemble.metNames      = metsList(2:end,2);
ensemble.rxnMets       = cell(length(ensemble.rxnNames),1);
ensemble.metsBalanced  = find(xMets(:,1));
ensemble.metsSimulated = find(xMets(:,2));
ensemble.metsFixed     = find(xMets(:,3));
ensemble.Sred          = ensemble.S(ensemble.metsBalanced,ensemble.activeRxns);   % Reduced stoichiometry for kinetic model simulation
ensemble.Sred(sum(abs(ensemble.Sred),2)==0,:) = [];                               % Remove zero rows
ensemble.Sred(sum(ensemble.Sred~=0,2)==1,:)   = [];                               % Remove unbalanced mets
disp('General information loaded.');

% Define exchanges and corresponding std
ensemble.measRates    = zeros(size(measRates,1),ensemble.numConditions);
ensemble.measRatesStd = ensemble.measRates;
for ix = 1:ensemble.numConditions+1                     % We have to add the reference state
    ensemble.measRates(:,ix)    = measRates(:,2*ix-1);
    ensemble.measRatesStd(:,ix) = measRates(:,2*ix);
end

% Add split ratios (if any)
ensemble.splitRatios = [];
if ~isempty(splitRatios)
    ensemble.splitRatios = splitRatios';
end

% Add pool constraints (if any)
ensemble.poolConst = [];
if ~isempty(poolConstraints)
    for ix = 1:size(poolConstraints,2)
        ensemble.poolConst{ix} = poolConstraints(:,ix)';
    end
end

% Add ineq. thermodynamic constraints (if any)
ensemble.ineqThermoConst = [];
if ~isempty(ineqConstraints)
    ensemble.ineqThermoConst = ineqConstraints';
end

%% 2. Define matrix for flux calculation and compute fluxes based on measured Rates
Sflux                 = ensemble.S(ensemble.metsBalanced,:);
ensemble.expFluxes    = zeros(length(ensemble.activeRxns),ensemble.numConditions);
ensemble.expFluxesStd = ensemble.expFluxes;
idxMeas               = idxMeas(2:end,1);                  % Extract only the meaningful information
for ix = 1:ensemble.numConditions+1
    xMean = zeros(size(Sflux,2),1);
    xStd  = xMean;
    for jx = 1:size(idxMeas)                               % Extract information in the right order
        index        = strcmp(ensemble.rxns,idxMeas(jx));
        xMean(index) = ensemble.measRates(jx,ix);
        xStd(index)  = ensemble.measRatesStd(jx,ix);
    end
    
    % Compute fluxes robustly (add split ratios if any)
    if robustFluxes
        if (~isempty(ensemble.splitRatios) && any(ensemble.splitRatios(ix,:)))
            [vMean,vStd] = computeRobustFluxes([Sflux;ensemble.splitRatios(ix,:)],xMean,xStd);
        else
            [vMean,vStd] = computeRobustFluxes(Sflux,xMean,xStd);
        end
    else
        vMean = xMean;
        vStd  = xStd;
    end
    
    % Calculate reference state and exp conditions fluxes
    if (ix==1)
%         dirFluxRef          = sign(vMean);                             % Define directionality at Ref for thermo calculations
        ensemble.fluxRef    = vMean(ensemble.activeRxns);              % Assign only the modelleded rxns
        ensemble.fluxRefStd = vStd(ensemble.activeRxns);
    else
        ensemble.expFluxes(:,ix-1)    = vMean(ensemble.activeRxns);
        ensemble.expFluxesStd(:,ix-1) = vStd(ensemble.activeRxns);
    end
end

% Get free fluxes
P = 1*(null(ensemble.Sred,'r')~=0);                                                             % Extract nullspace pattern, build adjacency matrix and pivots
ensemble.freeFluxes = [];
for ix = 1:size(P,2)
    for jx = 1:size(P,1)
        if (sum((P(jx,ix:end)~=0),2)==1)
            ensemble.freeFluxes = [ensemble.freeFluxes;jx];
            break;
        end
    end
end
ensemble.simWeights = ensemble.expFluxes(ensemble.freeFluxes,:);                                % Define simulation weights (these were based here on the flux magnitude)
disp('Flux data computed and loaded.');

%% 3. Perform thermodynamic calculations
if computeThermo
    idxNotExch  = find(~ismember(1:numel(ensemble.rxns),ensemble.exchRxns));
    Sthermo     = ensemble.S(:,idxNotExch);
    DGr_std     = xDG_std(idxNotExch,:);                                                        % Use only reactions with known thermodynamics
    vmin        = ensemble.fluxRef - 2*ensemble.fluxRefStd;
    vmax        = ensemble.fluxRef + 2*ensemble.fluxRefStd;
    xmin        = xMetsThermo(:,1);
    xmax        = xMetsThermo(:,2);
    DGr_std_min = DGr_std(:,1);
    DGr_std_max = DGr_std(:,2);
    gibbsRanges = computeGibbsFreeEnergyRanges(Sflux,Sthermo,DGr_std_min,DGr_std_max,vmin,vmax,xmin,xmax,idxNotExch,ineqConstraints', ensemble.rxns);     
    ensemble.gibbsRanges               = -1e2*ones(size(Sfull,1),2);                            % Allocate memory for DGr calculations
    ensemble.gibbsRanges(idxNotExch,:) = gibbsRanges;                                           % Remove thermodynamic info from exchang rxns
else
    ensemble.gibbsRanges = xDG_std;
end
disp('Thermodynamic data computed and loaded.');

%% 4. Load metabolomic data
% Determine the kinetically active rxns
idxMets = idxMets(2:end,1);               % Extract only the meaningful information
ensemble.metsDataMin  = zeros(length(ensemble.metsSimulated),ensemble.numConditions);
ensemble.metsDataMax  = ensemble.metsDataMin;
ensemble.metsDataMean = ensemble.metsDataMin;
for ix = 1:ensemble.numConditions
    for jx = 1:size(idxMets,1)            % Extract information in the right order
        index = strcmp(ensemble.mets(ensemble.metsSimulated),idxMets(jx));
        ensemble.metsDataMin(index,ix)  = metsData(jx,3*ix-2);
        ensemble.metsDataMean(index,ix) = metsData(jx,3*ix-1);
        ensemble.metsDataMax(index,ix)  = metsData(jx,3*ix);
    end
end
invalidMets                          = (ismember(ensemble.metsSimulated,ensemble.metsFixed));          % Keep only the information related to simulated mets not fixed
ensemble.metsDataMin(invalidMets,:)  = [];
ensemble.metsDataMean(invalidMets,:) = [];
ensemble.metsDataMax(invalidMets,:)  = [];
ensemble.metsDataMin(ensemble.metsDataMin==0)   = 1e-1;                                                % ~ 2 orders of magnitude of dynamic range
ensemble.metsDataMean(ensemble.metsDataMean==0) = 1e0;
ensemble.metsDataMax(ensemble.metsDataMax==0)   = 1e1;
disp('Metabolomics data loaded.');

%% 5. Extract kinetic information
ensemble.prevPrior = zeros(1,ensemble.numStruct);
ensemble.prevPriorInfo{ensemble.numStruct,2} = [];

for jx = 1:ensemble.numStruct
    try
        [priorRxns,priorParams] = xlsread(xlsxFile,['prior',num2str(jx)]);
        if ~isempty(priorRxns) && ~isempty(priorParams)
            ensemble.prevPrior(jx) = 1;
            ensemble.prevPriorInfo{jx,1} = priorRxns;
            ensemble.prevPriorInfo{jx,2} = priorParams;
        end
    catch
        disp(['No previous prior for model ',num2str(jx)]);
    end
    try
        [xKinetic,strKinetic] = xlsread(xlsxFile,['kinetics',num2str(jx)]);                    % read kinetic info from structure jx
        strKinetic = fixVariableNames(strKinetic, 'r', 'kinetics');     
    catch
        error("The kinetics sheet couln't be read. Make sure it is named as kinetics1.");
        break;
    end
    
    % Load kinetic information
    strKinetic(1,:) = [];                                                                      % remove useless information
    if (size(strKinetic,1)==length(ensemble.activeRxns))
        ensemble.allosteric{jx} = zeros(size(ensemble.activeRxns));
        ensemble.subunits{jx}   = zeros(size(ensemble.activeRxns));
        ensemble.rxnMechanisms{jx}{size(strKinetic,1),1} = [];
        ensemble.inhibitors{jx}{size(strKinetic,1),1}    = [];
        ensemble.activators{jx}{size(strKinetic,1),1}    = [];
        ensemble.negEffectors{jx}{size(strKinetic,1),1}  = [];
        ensemble.posEffectors{jx}{size(strKinetic,1),1}  = [];
        for ix = 1:size(strKinetic,1)
            index = find(ismember(ensemble.rxns(ensemble.activeRxns),strKinetic(ix)));
            ensemble.allosteric{jx}(index)    = xKinetic(ix,1);                                  % boolean vector indicating whether the enzyme is allosteric or not
            ensemble.subunits{jx}(index)      = xKinetic(ix,2);                                  % number of subunits present in the enzyme
            ensemble.rxnMechanisms{jx}{index} = strKinetic{ix,2};                                % string array with catalytic mechanisms for every reaction
            ensemble.subOrder{jx}{index}      = strKinetic{ix,3};                                % string array with the order of substrate bindings
            ensemble.prodOrder{jx}{index}     = strKinetic{ix,4};                                % string array with the order of product release
            ensemble.promiscuity{jx}{index}   = strKinetic{ix,5};                                % string array with rxn names this one is promiscuous with
            ensemble.inhibitors{jx}{index}    = strKinetic{ix,6};                                % string array with the inhibitor names for each reaction
            ensemble.activators{jx}{index}    = strKinetic{ix,7};                                % string array with the activator names for each reaction
            ensemble.negEffectors{jx}{index}  = strKinetic{ix,8};                                % string array with negative effector names for each reaction
            ensemble.posEffectors{jx}{index}  = strKinetic{ix,9};                                % string array with positive effector names for each reaction
        end
    else
        disp('The number of active rxns does not match the number of kinetic mechanisms.'); 
        break;
    end    
    if (jx==1)
        ensemble.kinActRxns   = find(~strcmp(ensemble.rxnMechanisms{jx},'fixedExchange'));
        
        %% 5/6. Extract proteomic information
        idxProt = idxProt(2:end,1);               % Extract only the meaningful information
        ensemble.protDataMin  = zeros(length(ensemble.kinActRxns),ensemble.numConditions);
        ensemble.protDataMax  = ensemble.protDataMin;
        ensemble.protDataMean = ensemble.protDataMin;
        for lx = 1:ensemble.numConditions
            for kx = 1:size(idxProt,1)            % Extract information in the right order
                index = strcmp(ensemble.rxns(ensemble.kinActRxns),idxProt(kx));
                ensemble.protDataMin(index,lx)  = protData(kx,3*lx-2);
                ensemble.protDataMean(index,lx) = protData(kx,3*lx-1);
                ensemble.protDataMax(index,lx)  = protData(kx,3*lx);
            end
        end
        rxnWithNoProteinData = find(all(ensemble.protDataMax'==1)&all(ensemble.protDataMin'==1)&all(ensemble.protDataMean'==1));        % These rxns are a different kind of 'inactive rxns' thus they are still considered in the active field
        ensemble.protDataMin(rxnWithNoProteinData,:)  = [];                                                                             % remove protein entries for rxns without protein information
        ensemble.protDataMax(rxnWithNoProteinData,:)  = [];
        ensemble.protDataMean(rxnWithNoProteinData,:) = [];
        disp('Proteomics and exchange data loaded and consistent.');
        
        % Figure out inactive rxns
        ensemble.kinInactRxns = find(strcmp(ensemble.rxnMechanisms{jx},'fixedExchange'));
        ensemble.fixedExch    = [ensemble.fluxRef(ensemble.kinInactRxns,:),ensemble.expFluxes(ensemble.kinInactRxns,:)];
        ensemble.kinInactRxns = [ensemble.kinInactRxns;rxnWithNoProteinData'];
        if ~isempty(ensemble.fixedExch) && ~isempty(rxnWithNoProteinData)
            ensemble.fixedExch = [ensemble.fixedExch;ones(numel(rxnWithNoProteinData),ensemble.numConditions+1)];
        elseif ~isempty(rxnWithNoProteinData)
            ensemble.fixedExch = ones(numel(rxnWithNoProteinData),ensemble.numConditions+1);
        end
    end          
        
    % Include information related to the activators and inhibitors, 
    % promiscuous reactions and substrate binding/release order
    for ix = 1:size(ensemble.rxnMechanisms{jx},1)
        if ~isempty(ensemble.promiscuity{jx}{ix})
           
            promiscuous_rxns_list  = regexp(ensemble.promiscuity{jx}{ix},' ','split');
            ensemble.promiscuity{jx}{ix} = find(ismember(ensemble.rxns(ensemble.activeRxns), promiscuous_rxns_list{1}));
            
            for rxn_i = 2:size(promiscuous_rxns_list, 2)
                ensemble.promiscuity{jx}{ix} = [ensemble.promiscuity{jx}{ix} find(ismember(ensemble.rxns(ensemble.activeRxns), promiscuous_rxns_list{rxn_i}))];
            end
            ensemble.promiscuity{jx}{ix} = sort(ensemble.promiscuity{jx}{ix});
        end
        if ~isempty(ensemble.subOrder{jx}{ix})
            ensemble.subOrder{jx}{ix}   = regexp(ensemble.subOrder{jx}{ix},' ','split');
        end
        if ~isempty(ensemble.prodOrder{jx}{ix})
            ensemble.prodOrder{jx}{ix}   = regexp(ensemble.prodOrder{jx}{ix},' ','split');
        end
        if ~isempty(ensemble.inhibitors{jx}{ix})
            ensemble.inhibitors{jx}{ix}   = regexp(ensemble.inhibitors{jx}{ix},' ','split');
        end
        if ~isempty(ensemble.activators{jx}{ix})
            ensemble.activators{jx}{ix}   = regexp(ensemble.activators{jx}{ix},' ','split');
        end
        if ~isempty(ensemble.negEffectors{jx}{ix})
            ensemble.negEffectors{jx}{ix} = regexp(ensemble.negEffectors{jx}{ix},' ','split');
        end
        if ~isempty(ensemble.posEffectors{jx}{ix})
            ensemble.posEffectors{jx}{ix} = regexp(ensemble.posEffectors{jx}{ix},' ','split');
        end
    end
    
    % Read each mechanism and write the reaction patterns
    ensemble.kineticFxn{jx} = [ensemble.description,'_Kinetics',num2str(jx)];      % string with kinetic function name
    mkdir('reactions');
    for ix = 1:size(ensemble.rxnMechanisms{jx},1)
        
        % Non-enzymatic reactions: fixed exchange
        if strcmp(ensemble.rxnMechanisms{jx}{ix},'fixedExchange')
            buildFixedExchange(ensemble.rxns{ix},jx)
            ensemble.revMatrix{ix,jx}   = [];
            ensemble.forwardFlux{ix,jx} = [];
            ensemble.Nelem{ix,jx}       = [];
        
        % Non-enzymatic reactions: free exchange
        elseif strcmp(ensemble.rxnMechanisms{jx}{ix},'freeExchange')
            buildFreeExchange(ensemble.rxns{ix},jx)
            ensemble.revMatrix{ix,jx}   = [];
            ensemble.forwardFlux{ix,jx} = [];
            ensemble.Nelem{ix,jx}       = [];            
            
        % Non-enzymatic reactions: passive difusion
        elseif strcmp(ensemble.rxnMechanisms{jx}{ix},'diffusion')
            buildDiffusion(ensemble.rxns{ix},jx)
            ensemble.revMatrix{ix,jx}   = [];
            ensemble.forwardFlux{ix,jx} = [];
            ensemble.Nelem{ix,jx}       = [];
        
        % Non-enzymatic reactions: mass action
        elseif strcmp(ensemble.rxnMechanisms{jx}{ix},'massAction')
            buildMassAction(ensemble.rxns{ix},jx)
            ensemble.revMatrix{ix,jx}   = [];
            ensemble.forwardFlux{ix,jx} = [];
            ensemble.Nelem{ix,jx}       = [];
            
            % Enzymatic reactions
        else
            
            if size(ensemble.promiscuity{jx}{ix}) > 0
                promiscuousRxnI = find(ensemble.promiscuity{jx}{ix} == ix);
            else 
               promiscuousRxnI = 0; 
            end

            % Allosteric enzymes
            if ensemble.allosteric{jx}(ix)
                [revMatrix,forwardFlux,metList] = reactionPattern(ensemble.rxnMechanisms{jx}{ix},ensemble.rxns{ix},2,jx, promiscuousRxnI);
                buildAllosteric(metList,[ensemble.rxns{ix},num2str(jx)],ensemble.negEffectors{jx}{ix},ensemble.posEffectors{jx}{ix})
                ensemble.metLists{ix,1} = metList;
                % Non-allosteric enzymes
            else
                [revMatrix,forwardFlux, metList] = reactionPattern(ensemble.rxnMechanisms{jx}{ix},ensemble.rxns{ix},1,jx, promiscuousRxnI);
                ensemble.metLists{ix,1} = metList;
            end
                        
            % Build Selem based on the mechanism stoichiometry
            Selem = [];
            for i = 1:max(forwardFlux(:))
                Selem = [Selem;((forwardFlux(:,1)==i)-(forwardFlux(:,2)==i))'];
            end
            
            % Save features of the kinetic mechanism
            ensemble.revMatrix{ix,jx}   = revMatrix;               % reversibility matrix for the reaction mechanism
            ensemble.forwardFlux{ix,jx} = forwardFlux;             % forward flux for the reaction mechanism
            
            % Compute null(Selem) for branching calculations and precondition accordingly
            Ntemp    = null(Selem,'r');
            ixNnzCol = find(any(Ntemp<0),1,'first');
            
            % Check that Nelem is suitable for branching calculations (loop until columns have only nonnegative entries)
            while ~isempty(ixNnzCol)
                for cx = 1:size(Ntemp,2)
                    if (cx==ixNnzCol)
                        continue;
                    else
                        Npos = Ntemp(:,ixNnzCol)+Ntemp(:,cx);
                        Nneg = Ntemp(:,ixNnzCol)-Ntemp(:,cx);
                        if ~any(Npos<0)
                            Ntemp(:,ixNnzCol) = Npos; break;
                        elseif ~any(Nneg<0)
                            Ntemp(:,ixNnzCol) = Nneg; break;
                        end
                    end
                end
                ixNnzCol = find(any(Ntemp<0),1,'first');
            end
            ensemble.Nelem{ix,jx} = Ntemp;                         % assign final positive null basis
        end
    end
    copyfile('./reactions',['./reactions',num2str(jx)]);
    addpath(['./reactions',num2str(jx)]);
    rmdir('./reactions','s');
    
    % Build kinetic fxn and find active species (do not build hess partern)
    [rxnMetLinks,freeVars,metsActive] = buildKineticFxn(ensemble,ensemble.kineticFxn{jx},jx);
    if (jx==1)                                                 % freeVars and freeMets are the same for all the structure
        ensemble.rxnMetLinks = rxnMetLinks;
        ensemble.freeVars    = freeVars;
        ensemble.metsActive  = metsActive;
    end    
    buildKineticFxn(ensemble,ensemble.kineticFxn{jx},jx);
    disp(['Kinetic information loaded and kinetic model built: Structure ',num2str(jx),'.']);
end

disp('Ensemble structure ready.');


