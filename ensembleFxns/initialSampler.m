function [models,strucIdx,xopt,tolScore,simulatedFlux] = initialSampler(ensemble)
%--------------------------------------------------------------------------
% Sampling initial ensemble of kinetic models using rejection scheme
%
% Inputs:       ensemble (structure) , workerIdx (double)
%
% Outputs:      initial sampled ensemble structure
%--------------------- Pedro Saa 2016 -------------------------------------
%% Initialze parameters
RT       = 8.314*298.15/1e3;                                               % gas constant times the absolute temperature (298.15 K)
massTol  = size(ensemble.Sred,1)*1e-10;								       % #balances*tol^2

% Figure out NLP solver
if strcmpi(ensemble.solver,'NLOPT')											   % Solver parameters for NLOPT
    opt.algorithm = 40; 									   			   % 11(NLOPT_LD_LBFGS), options: 40(NLOPT_LD_SLSQP), 13(NLOPT_LD_VAR1), 14(NLOPT_LD_VAR2)
    opt.ftol_abs  = 1e-11;
    opt.xtol_abs  = 1e-10*ones(1,numel(ensemble.freeVars));
    opt.maxeval   = 1e4;
elseif strcmpi(ensemble.solver,'FMINCON') || isempty(ensemble.solver)							   % Solver parameters (FMINCON)
    options = optimset('Display','off','Algorithm','sqp','MaxIter',1e4,'TolFun',1e-11,'TolX',1e-10);
end

% Check if there are pool constraints
if ~isempty(ensemble.poolConst)
    for ix = 1:numel(ensemble.poolConst)
        A{ix} = ensemble.poolConst{ix}(1:numel(ensemble.metsActive));      % extract rhs of from pool constraint matrix
        b{ix} = ensemble.poolConst{ix}(numel(ensemble.metsActive)+1:end);
    end
else
    A = [];                                                                % inequality constraints matrix
    b = [];                                                                % rhs ineequality constraints
end
Aeq    = [];                                                               % equality constraints matrix
beq    = [];                                                               % rhs equality constraints
x0     = [ensemble.metsDataMean;ensemble.protDataMean];                    % initial guess
lb     = [ensemble.metsDataMin;ensemble.protDataMin];                      % lower bounds on free vars
ub     = [ensemble.metsDataMax;ensemble.protDataMax];                      % upper bounds on free vars
nlcons = [];                                                               % nonlinear constraints (not used)

%% Execute Rejection-ABC
acceptanceRate = 1;
counter        = 0;

% Loop until the number of valid particles is reached
while true

    % Update attempt counter
    counter = counter+1;

    % Sample model structure
    strucIdx = randi(ensemble.numStruct);

    % Sample pool parameters (if any)
    if ~isempty(ensemble.poolConst)
        poolFactor{numel(ensemble.poolConst)} = [];
        for ix = 1:numel(ensemble.poolConst)

            % Generate pool factor ~ Dir(alpha) using independent gamma distributions
            alphaPoolFactor = ensemble.populations(1).probParams(strucIdx).alphaPoolFactor{ix};
            poolFactorTemp  = randg(alphaPoolFactor);
            poolFactorTemp  = poolFactorTemp/sum(poolFactorTemp);

            % Update pool constraint matrix accordingly
            A_opt{ix} = A{ix};
            A_opt{ix}(A{ix}~=0) = poolFactorTemp;

            % Save sampled poolfactor
            poolFactor{ix} = poolFactorTemp;
        end
        models(1).poolFactor = poolFactor;
    else
        models(1).poolFactor = [];
    end

    % Randomly distribute flux between isoenzymes

    if ~isempty(ensemble.uniqueIso)
        for xi = 1:size(ensemble.uniqueIso,1)
            group = find(strcmp(ensemble.isoenzymes,ensemble.uniqueIso{xi}));
            splitFactor = zeros(size(group,1),1);
            totalFlux = sum(ensemble.fluxRef(group));
            for yi = 1:size(splitFactor,1)
                splitFactor(yi) = randg();
            end
            splitFactor = splitFactor./sum(splitFactor);
            ensemble.fluxRef(group) = splitFactor.*totalFlux;
        end
    end

    models.refFlux =  ensemble.fluxRef;

    % Sample gibbs free energy of reactions
    gibbsFactor = mvnrnd(ensemble.populations(1).probParams(strucIdx).muGibbsFactor,ensemble.populations(1).probParams(strucIdx).sigmaGibbsFactor)';
    gibbsFactor = exp(gibbsFactor)./(1 + exp(gibbsFactor));
    gibbsEnergy = gibbsFactor.*ensemble.gibbsRanges(ensemble.thermoActive,2) + (1-gibbsFactor).*ensemble.gibbsRanges(ensemble.thermoActive,1);
    models(1).gibbsFactor = gibbsFactor;

    % Determine gibbs free energy of reaction
    ensemble = sampleGibbsReactionEnergies(ensemble,gibbsEnergy,strucIdx);

    % Sample Reversibilities
    [ensemble, models] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx);

    % Sample enzyme abundances
    models = sampleEnzymeAbundances(ensemble,models,strucIdx);

    % Sample modifier elementary fluxes (positions are given where exp(R)=1)
    models = sampleModifierElemFluxes(ensemble, models, strucIdx);

    % Calculate rate parameters for allosteric reaction part;
    [ensemble, models] = sampleAllostery(ensemble, models, strucIdx);

    %thermoCounter   = 1;
    for activRxnIdx = 1:numel(ensemble.kinActRxns)
        %disp(ensemble.rxns(ensemble.kinActRxns(activRxnIdx)));

        % Case 1: Diffusion and Exchanges
        if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
                strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')
            models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).kineticParams = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));

            % Case 2: Enzymatic reactions
        else

            % Check whether the reaction is mass action
            if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
                reactionFlux = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));
                gibbsTemp =  ensemble.gibbsTemp{ensemble.kinActRxns(activRxnIdx)};
                models(1).rxnParams(activRxnIdx).kineticParams = [1,exp(gibbsTemp/RT)]*reactionFlux/(1-exp(gibbsTemp/RT));
                continue;
            end

            promiscRxnsList = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
            revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
            reverTemp = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};
            reactionFlux = ensemble.reactionFluxAllosteric(ensemble.kinActRxns(activRxnIdx));
            randomEnzymesR = models(1).rxnParams(activRxnIdx).enzymeAbundances';
            extremePathways = ensemble.extremePathways{strucIdx}{activRxnIdx};


            % Sample branching factor (if necessary)
            branchFactor = 1;
            Nelem        = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};

            % If the reaction is promiscuous
            if size(promiscRxnsList,1) > 0
                branchFactor = ensemble.reactionFluxAllosteric(promiscRxnsList);

                % If the promiscuous reactions share common steps
                if sum(sum(Nelem)) > size(Nelem,1)
                    reactionFlux = sum(ensemble.reactionFluxAllosteric(promiscRxnsList));

                % If the promiscuous reactions do not share any common steps
                else
                    reactionFlux = max(ensemble.reactionFluxAllosteric(promiscRxnsList));
                end

                % For non promiscuous reactions
            else
                if (size(extremePathways,2)>1)
                    branchFactor = zeros(1,size(extremePathways,2));
                    for ix = 1:size(extremePathways,2)
                        aBranch            = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:));
                        branchFactor(1, ix) = aBranch;
                    end
                    branchFactor = branchFactor/sum(branchFactor);
                end
            end
            models(1).rxnParams(activRxnIdx).branchFactor = branchFactor';

            % Get modifier elementary fluxes (positions are given were exp(R)=1)
            modifierElemFlux = models(1).rxnParams(activRxnIdx).modiferElemFlux';

            % VI. Calculate rate parameters
            forwardFlux    = ensemble.forwardFlux{ensemble.kinActRxns(activRxnIdx),strucIdx};
            models(1).rxnParams(activRxnIdx).kineticParams = ...
                calculateKineticParams(reverTemp,forwardFlux,reactionFlux,randomEnzymesR,extremePathways,branchFactor,modifierElemFlux);
        end
    end

    % Test model consistency
    kineticFxn = str2func(ensemble.kineticFxn{strucIdx});
    testFlux   = feval(kineticFxn,ones(size(ensemble.freeVars,1),1),models,ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);

    % If the model is consistent continue
    condition = true;
    if any(abs(testFlux-ensemble.fluxRef)>1e-6)
        condition = false;
        disp(['There are consistency problems during the reaction sampling. Model ID: ',num2str(strucIdx)]);
    end

    % Check sampling mode. For the ORACLE mode, no need to simulate
    if strcmpi(ensemble.sampler,'ORACLE'); break;

        % For the remaining modes, we need to simulate the model in the
        % experimental conditions
    elseif ~strcmpi(ensemble.sampler,'ORACLE') && condition

        % Simulate fluxes
        tolScore      = [];
        simulatedFlux = zeros(numel(ensemble.activeRxns),ensemble.numConditions);
        xopt          = zeros(size(x0,1),ensemble.numConditions);

        % Solver call (OPTI Toolbox not implemented yet)
        for ix = 1:ensemble.numConditions

            % NLOPT call
            if strcmpi(ensemble.solver,'NLOPT')

                % Define anonymous function with objective and bound
                % constraints
                opt.min_objective = @(x) kineticFxn(x,models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
                opt.lower_bounds  = lb(:,ix);
                opt.upper_bounds  = ub(:,ix);

                % Solve S*v(k,X) = 0; s.t. A*X <= beq, lb < X <ub, with extra constraints (e.g., pool or ratio constraints). Otherwise solve solve S*v(k,X) = 0; s.t. lb < X <ub, with no extra constraints
                if ~isempty(ensemble.poolConst)
                    for jx = 1:numel(ensemble.poolConst)
                        opt.fc{1,2*jx-1} = (@(x) poolConstraintFxn(x,[A_opt{jx},zeros(1,numel(ensemble.freeVars)-numel(ensemble.metsActive))],b{jx}(2*ix-1)));
                        opt.fc{1,2*jx}   = (@(x) poolConstraintFxn(x,[-A_opt{jx},zeros(1,numel(ensemble.freeVars)-numel(ensemble.metsActive))],-b{jx}(2*ix)));
                    end
                    opt.fc_tol = 1e-6*ones(1,2*numel(ensemble.poolConst));
                end
                [xopt(:,ix),fmin] = nlopt_optimize(opt,x0(:,ix));

                % FMINCON call
            else
                [xopt(:,ix),fmin] = fmincon(kineticFxn,x0(:,ix),[],[],[],[],lb(:,ix),ub(:,ix),[],options,models,ensemble.fixedExch,ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);

                % Pool constraints not implemented yet for this solver
            end

            % Check mass balance consistency
            if (fmin<massTol)

                % Simulate fluxes if the system is mass-balanced
                simulatedFlux(:,ix) = feval(kineticFxn,xopt(:,ix),models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);

                % Calculate discrepancy score
                tolScore = [tolScore,max(sqrt(mean(((simulatedFlux(ensemble.freeFluxes,ix)-ensemble.simWeights(:,ix))./ensemble.simWeights(:,ix)).^2)))];

                % Check tolerance inmediately for this condition
                if (tolScore(end)>ensemble.tolerance(1))
                    condition = false;
                    break;
                end
            else
                condition = false;
                break;
            end
        end

        % Compute tolerance, acceptance rate and break
        if condition
            tolScore       = max(tolScore);                 % Infinity norm as discrepancy measure
            acceptanceRate = 1/counter;
            break;

            % Delete model if not accurate or mass-balanced
        else
            models = [];
        end
    end
end

% Save results and write progress to a temp file (except for the ORACLE mode)
if ~strcmpi(ensemble.sampler,'ORACLE')
    try
        load progress.txt
        progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
        save progress.txt -ascii progress
        save(['temp/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux');

        % If another worker is writing on the file, wait a brief random time
    catch
        pause(randi(2)*rand(1));
        load progress.txt
        progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
        save progress.txt -ascii progress
        save(['temp/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux');
    end
end
