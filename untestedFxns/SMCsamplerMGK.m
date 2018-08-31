function [models,strucIdx,xopt,tolScore,simulatedFlux,log_particleWeight] = SMCsamplerMGK(ensemble,popIdx,currTol,log_weightThreshold,prevParams)
%--------------------------------------------------------------------------
% Sampling ensemble of kinetic parameters using adaptive SMC
%
% Inputs:       ensemble (structure)
%
% Outputs:      initial sampled ensemble structure
%--------------------- Pedro Saa 2016 -------------------------------------
%% Initialze parameters
RT       = 8.314*298.15/1e3;                                               % gas constant times the absolute temperature (298.15 K)
massTol  = size(ensemble.Sred,1)*1e-9;								       % #balances*tol^2
DDG_max  = 0;															   % max. Gibbs free energy difference of conformation (kJ/mol)
DDG_min  = -25;															   % min. Gibbs free energy difference of conformation (kJ/mol) ~ 1 ATP hydrolysis
bind_max = 50;															   % max. fraction of Keff/Eff
bind_min = .1;															   % min. fraction of Keff/Eff
log_maxWeightDiff = 10;													   % sampling parameter (optional), controls the maximum importance weight

% Optimizer parameters
opt.algorithm = 40; 										   			   % 11(NLOPT_LD_LBFGS), options: 40(NLOPT_LD_SLSQP), 13(NLOPT_LD_VAR1), 14(NLOPT_LD_VAR2)
opt.ftol_abs  = 1e-11;
opt.xtol_abs  = 1e-10*ones(1,numel(ensemble.freeVars));
opt.maxeval   = 5e3;
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
beq    = [];
x0     = [ensemble.metsDataMean;ensemble.protDataMean];                    % initial guess
lb     = [ensemble.metsDataMin;ensemble.protDataMin];                      % lower bounds on free vars
ub     = [ensemble.metsDataMax;ensemble.protDataMax];                      % upper bounds on free vars
nlcons = [];                                                               % nonlinear constraints (not used)

%% Execute Adaptive ABC-SMC sampler
acceptanceRate = 1;
counter        = 0;

% Loop until the number of valid particles is reached
while true
    
    % Update attempt counter
    counter = counter+1;
    
    % Sample model structure indicator
    strucIdx = sum(rand(1) >= cumsum([0,ensemble.populations(end-1).structureWeights]));
    
    % Sample an alive particle from the previous population based on the structure weights
    activeWeights = ensemble.populations(popIdx-1).weights(ensemble.populations(popIdx-1).strucIdx==strucIdx);
    paramIdx      = sum(rand(1) >= cumsum([0;activeWeights]/ensemble.populations(end-1).structureWeights(strucIdx)));
    
    % Perturb the chosen particle
    sigma         = ensemble.populations(popIdx).probParams(strucIdx).globalKernelParam;									% Extract model-specific cov matrix
    proposalArray = mvnrnd(prevParams(strucIdx).paramArray(:,paramIdx)',sigma);												% NOTE: the rejection could be done inmediately if wanted, however the proposed array must be converted to a model structure
    indexes       = prevParams(strucIdx).indexes;
    
	% Get index of the previous model
	activeParticles = find(ensemble.populations(popIdx-1).strucIdx==strucIdx);
	modelIdx        = activeParticles(paramIdx);
	
    % Initialize parameter coordinate counter
    coordIdx = 1;
    
    % Sample pool parameters (if any)
    if ~isempty(ensemble.poolConst)
        poolFactor{numel(ensemble.poolConst)} = [];
        for ix = 1:numel(ensemble.poolConst)
            
            % Map back to the original space
            poolFactorTemp = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
            poolFactorTemp = map2simplex(poolFactorTemp);
            
            % Update pool constraint matrix accordingly
            A_opt{ix} = A{ix};
            A_opt{ix}(A{ix}~=0) = poolFactorTemp;
                        
            % Save sampled poolfactor
            poolFactor{ix} = poolFactorTemp;
            coordIdx = coordIdx + 1;                                                                                                                % advance parameter counter
        end
        models(1).poolFactor = poolFactor;
    else
        models(1).poolFactor = [];
    end
    
    % Form gibbs free energy vector
    gibbsFactor = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
    gibbsFactor = (exp(gibbsFactor)./(1 + exp(gibbsFactor)))';																						% apply inverse logit transformation
    gibbsEnergy = gibbsFactor.*ensemble.gibbsRanges(ensemble.thermoActive,2) + (1-gibbsFactor).*ensemble.gibbsRanges(ensemble.thermoActive,1);
    coordIdx    = coordIdx + 1;																														% advance parameter counter
    models(1).gibbsFactor = gibbsFactor';
    
    % Calculation of reversibilities for each raction in the model
    thermoCounter   = 1;
    for activRxnIdx = 1:numel(ensemble.kinActRxns)
        
        % Case 1: Diffusion and Exchanges
        if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'diffusion')||...
                strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'freeExchange')
            models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).kineticParams = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));
            
            % Case 2: Enzymatic reactions
        else
            
            % Determine gibbs free energy of reaction
            if ismember(ensemble.kinActRxns(activRxnIdx),ensemble.thermoActive)
                gibbsTemp     = gibbsEnergy(thermoCounter);
                thermoCounter = thermoCounter + 1;
            elseif (ensemble.gibbsRanges(activRxnIdx,1)==ensemble.gibbsRanges(activRxnIdx,2))
                gibbsTemp     = ensemble.gibbsRanges(activRxnIdx,1);
            end
            
            % Check whther the reaction is mass action
            if strcmp(ensemble.rxnMechanisms{strucIdx}{activRxnIdx},'massAction')
                reactionFlux = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));
                models(1).rxnParams(activRxnIdx).kineticParams = [1,exp(gibbsTemp/RT)]*reactionFlux/(1-exp(gibbsTemp/RT));
                continue;
            end
            
            % Form new reversibility (Note: this is the normalized parameter)
            revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
            randomRev = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
            randomRev = map2simplex(randomRev)';																			% map onto the simplex
            if (size(revMatrix,1)>1)                                                                                     	% Determine whether there are branches for the reversibility calculation
                lbRev     = randomRev;
                randomRev = computeBranchedReversibilities(revMatrix,lbRev);
                models(1).rxnParams(activRxnIdx).reversibilities = lbRev';        										 	% Save transpose (note this has to be the lower bound!)
            else
                if (size(revMatrix,1)==1)&&(any(revMatrix==0))																	% Fix reversibilities with elementary fluxes
                    randomRev(revMatrix~=0) = randomRev;
                    randomRev(revMatrix==0) = 0;
                end
                models(1).rxnParams(activRxnIdx).reversibilities = randomRev';
            end
            coordIdx  = coordIdx + 1;																		 				% advance coordinate counter
            reverTemp = exp(randomRev*gibbsTemp/RT);                                                                     	% Convert to the proper units for later calculation
            
            % Form new enzyme abundances vector
            forwardFlux    = ensemble.forwardFlux{ensemble.kinActRxns(activRxnIdx),strucIdx};
            randomEnzymesR = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
            randomEnzymesR = map2simplex(randomEnzymesR)';																 	% map onto the simplex
            coordIdx       = coordIdx + 1;
            models(1).rxnParams(activRxnIdx).enzymeAbundances = randomEnzymesR';
            
            % II. Sample branching factor (if necessary)
            branchFactor = 1;
            Nelem        = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
            if (size(Nelem,2)>1)
                branchFactor = zeros(size(Nelem,2),1);
                for ix = 1:size(Nelem,2)
                    tempFactor       = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
                    branchFactor(ix) = exp(tempFactor)./(exp(tempFactor) + 1);													% apply inverse logit transformation
                    coordIdx   = coordIdx + 1;
                end
            end
            models(1).rxnParams(activRxnIdx).branchFactor = branchFactor';
            
            % III. Sample modifier elementary fluxes
            modifierElemFlux = [];
            if (size(revMatrix,1)==1)&&any(revMatrix==0)
                modifierElemFlux = zeros(sum(reverTemp==1),1);
                for ix = 1:sum(reverTemp==1)
                    tempFactor = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
                    modifierElemFlux(ix) = exp(tempFactor)./(exp(tempFactor) + 1);
                    coordIdx   = coordIdx + 1;
                end
            end
            models(1).rxnParams(activRxnIdx).modiferElemFlux = modifierElemFlux';      % save transpose of mod elem flux
            
            % IV. Calculate rate parameters
            reactionFlux = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));
            
            % If the reaction is allosteric
            if ensemble.allosteric{strucIdx}(ensemble.kinActRxns(activRxnIdx))
                
                % I. Sample allosteric parameters
                allostericParams = proposalArray(indexes(coordIdx,1):indexes(coordIdx,2));
                allostericParams = exp(allostericParams)./(1 + exp(allostericParams));									    	% apply inverse lgit transformation
                models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).allostericFactors = allostericParams;                 	% allosteric factors used in the calculation of L, KposEff & KnegEff
                subunits = ensemble.subunits{strucIdx}(ensemble.kinActRxns(activRxnIdx));
                coordIdx = coordIdx + 1;
                
                % II. Back-calculate kinetic parameters (order: L, posEff & negEff) (in the case of DXS proceed differently)
                if (ensemble.kinActRxns(activRxnIdx)~=1)
                    DDG = DDG_min + allostericParams(1)*(DDG_max - DDG_min);
                    L   = exp(-DDG/RT);
                    models(1).rxnParams(activRxnIdx).L = L;                                                                     % Save allosteric constant
                    
                    % III. Compute regulatoy term
                    Q = L*(randomEnzymesR(1))^subunits;
                    if (numel(allostericParams)>1)																				% There are possitive and/or negative effectors
                        posEffectors = ensemble.posEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)};                       % Possitive effectors
                        negEffectors = ensemble.negEffectors{strucIdx}{ensemble.kinActRxns(activRxnIdx)};                       % Negative effectors
                        
                        % IV. Back-calculate effector parameters
                        Keff = bind_min*exp(allostericParams(2:end)*log(bind_max/bind_min));
                        
                        % Both effectors present
                        if ~isempty(posEffectors) && ~isempty(negEffectors)
                            KposEff = Keff(1:max(size(posEffectors)));
                            KnegEff = Keff(max(size(posEffectors))+1:max(size(posEffectors))+max(size(negEffectors)));
                            models(1).rxnParams(activRxnIdx).KposEff = KposEff;							       		            % transpose params
                            models(1).rxnParams(activRxnIdx).KnegEff = KnegEff;
                            Q = Q*((1 + sum(1./KnegEff))/((1 + sum(1./KposEff))))^subunits;
                        elseif ~isempty(posEffectors)
                            models(1).rxnParams(activRxnIdx).KposEff = Keff;							                        % transpose params
                            Q = Q*((1 + sum(1./Keff)))^-subunits;
                        elseif ~isempty(negEffectors)
                            models(1).rxnParams(activRxnIdx).KnegEff = Keff;                        							% transpose params
                            Q = Q*((1 + sum(1./Keff)))^subunits;
                        end
                    end
                    regContribution = 1/(1 + Q);
                    
                    % For the special DXS case
                else
                    Ktpp  = (10*bind_min)*exp(allostericParams(1)*log((bind_max/5)/(10*bind_min)));								% modified allosteric constant
                    L     = (1 + 1/Ktpp)^-subunits;
                    Keff  = bind_min*exp(allostericParams(2:end)*log(bind_max/bind_min));
                    models(1).rxnParams(activRxnIdx).L = L;
                    models(1).rxnParams(activRxnIdx).KnegEff = Keff;
                    Q     = L*sum(1./Keff)^subunits;
                    regContribution = (1/(1 + Q))/subunits;
                end
                
                % V. Calculate regulatory and catalytic contributions at the reference state
                catContribution = reactionFlux/regContribution;
                reactionFlux    = catContribution/subunits;
            end
            
            % VI. Calculate rate parameters
            models(1).rxnParams(activRxnIdx).kineticParams = ...
                calculateKineticParams(reverTemp,forwardFlux,reactionFlux,randomEnzymesR,Nelem,branchFactor,modifierElemFlux);
        end
    end
    
    % Test model consistency
    kineticFxn = str2func(ensemble.kineticFxn{strucIdx});
    testFlux   = feval(kineticFxn,ones(size(ensemble.freeVars,1),1),models,ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);
    
    % If the model is consistent continue
    if all(abs(testFlux-ensemble.fluxRef)<1e-10)

        % Compute un-normalized particle weight and simulate only if greater than threshold (PRC step)
        log_particleWeight = computeParticleWeightsVectorizedMGK(models(1),strucIdx,ensemble,popIdx,prevParams);
        log_paccept        = min([0,log_particleWeight-log_weightThreshold]);
        
        if (rand(1)<exp(log_paccept)) && (log_particleWeight-log_weightThreshold<log_maxWeightDiff)

            % Update propsal weight
            log_particleWeight = log_particleWeight - log_paccept;
            
            % Simulate fluxes
            tolScore      = [];
            condition     = true;
            simulatedFlux = zeros(length(ensemble.activeRxns),ensemble.numConditions);
            for ix = 1:ensemble.numConditions
                
                % Define anonymous function
                opt.min_objective = @(x) kineticFxn(x,models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
                opt.lower_bounds  = lb(:,ix);
                opt.upper_bounds  = ub(:,ix);                 
                 
                % Solve S*v(k,X) = 0; s.t. A*X <= beq, lb < X <ub, with extra constraints (e.g., pool or ratio constraints). Otherwise solve solve S*v(k,X) = 0; s.t. lb < X <ub, with no extra constraints
                if ~isempty(ensemble.poolConst)
                    for jx = 1:numel(ensemble.poolConst)					
						opt.fc{1,2*jx-1} = (@(x) poolConstraintFxn(x,[A_opt{jx},zeros(1,numel(ensemble.freeVars)-numel(ensemble.metsActive))],b{jx}(2*ix-1)));
						opt.fc{1,2*jx}   = (@(x) poolConstraintFxn(x,[-A_opt{jx},zeros(1,numel(ensemble.freeVars)-numel(ensemble.metsActive))],-b{jx}(2*ix)));					
					end                
					opt.fc_tol  = 1e-6*ones(1,2*numel(ensemble.poolConst));
                end
				[xopt(:,ix),fmin]   = nlopt_optimize(opt,ensemble.populations(popIdx-1).xopt{modelIdx}(:,ix));
                simulatedFlux(:,ix) = feval(kineticFxn,xopt(:,ix),models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);
                
                % Check mass balance consistency
                if (fmin<massTol)
                    tolScore = [tolScore,max(sqrt(mean(((simulatedFlux(ensemble.freeFluxes,ix)-ensemble.simWeights(:,ix))./ensemble.simWeights(:,ix)).^2)))];
                    
                    % Check tolerance inmediately for this condition
                    if (tolScore(end)>currTol)
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
                tolScore       = max(tolScore);
                acceptanceRate = 1/counter;
                break;
                
                % Erase generated structure if it is not mass balanced
            else
                models = [];
            end
        else
            models = [];
        end
    end
end

% Save results and write progress to file
try
	load progress.txt
	progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
	save progress.txt -ascii progress
	save(['tempSMC/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux','log_particleWeight');

% If another worker is writing on the file, wait a small random time
catch
	pause(randi(2)*rand(1));
	load progress.txt
	progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
	save progress.txt -ascii progress
	save(['tempSMC/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux','log_particleWeight');
end