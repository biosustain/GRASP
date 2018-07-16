function [models,strucIdx,xopt,tolScore,simulatedFlux] = rejectionSampler(ensemble)
%--------------------------------------------------------------------------
% Sampling initial ensemble of kinetic models using rejection scheme
%
% Inputs:       ensemble (structure) , workerIdx (double)
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

% % Solver parameters (NLOPT)
% opt.algorithm = 40; 										   			   % 11(NLOPT_LD_LBFGS), options: 40(NLOPT_LD_SLSQP), 13(NLOPT_LD_VAR1), 14(NLOPT_LD_VAR2)
% opt.ftol_abs  = 1e-11;
% opt.xtol_abs  = 1e-10*ones(1,numel(ensemble.freeVars));
% opt.maxeval   = 1e4;

% Solver parameters (FMINCON)
options = optimset('Display','off','Algorithm','sqp','MaxIter',1e4,'TolFun',1e-11,'TolX',1e-10);
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
	
	% Sample gibbs free energy of reactions
	gibbsFactor = mvnrnd(ensemble.populations(1).probParams(strucIdx).muGibbsFactor,ensemble.populations(1).probParams(strucIdx).sigmaGibbsFactor)';
	gibbsFactor = exp(gibbsFactor)./(1 + exp(gibbsFactor));	
	gibbsEnergy = gibbsFactor.*ensemble.gibbsRanges(ensemble.thermoActive,2) + (1-gibbsFactor).*ensemble.gibbsRanges(ensemble.thermoActive,1);
	models(1).gibbsFactor = gibbsFactor;
	
    
    % Determine gibbs free energy of reaction
    ensemble = sampleGibbsReactionEnergies(ensemble, gibbsEnergy, strucIdx);
    
    % Sample Reversibilities
    [ensemble, models] = sampleGeneralReversibilities(ensemble, models, RT, strucIdx);
    
    % Calculate rate parameters for allosteric reaction part;
    [ensemble, models] = sampleAllostery(ensemble, models, strucIdx);
    
    % Sample modifier elementary fluxes (positions are given were exp(R)=1)
    [models] = sampleModifierElemFluxes(ensemble, models, strucIdx);
    
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
            
            promisc_rxns_list = ensemble.promiscuity{strucIdx}{ensemble.kinActRxns(activRxnIdx)};
            revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
            reverTemp = ensemble.reverTemp{ensemble.kinActRxns(activRxnIdx)};
            reactionFlux = ensemble.reactionFluxAllosteric(ensemble.kinActRxns(activRxnIdx));
            
            
            % B. Sample enzyme abundances
            alphaEnzymesR  = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances;
            
            % If it's a promiscuous reaction and it' not the first one in the set of promiscuous reactions
            if size(promisc_rxns_list,1) > 0 && ensemble.kinActRxns(activRxnIdx) ~= promisc_rxns_list(1)                    % Get the abundances from the first reaction in the set
                randomEnzymesR = models(1).rxnParams(promisc_rxns_list(1)).enzymeAbundances';
            else
                randomEnzymesR = randg(alphaEnzymesR');                                                                      % Sample from the gamma distribution with parameter alpha
                randomEnzymesR = randomEnzymesR/sum(randomEnzymesR);
            end
            models(1).rxnParams(activRxnIdx).enzymeAbundances = randomEnzymesR';
      
            % C. Sample branching factor (if necessary)
            branchFactor = 1;
            Nelem        = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};      
            
            % If the reaction is promiscuous 
            if size(promisc_rxns_list,1) > 0
                branchFactor = ensemble.reactionFluxAllosteric(promisc_rxns_list);
                if sum(sum(Nelem)) > size(Nelem,1)
                    reactionFlux = sum(ensemble.reactionFluxAllosteric(promisc_rxns_list));   % promiscuous reactions share common steps
                else
                    reactionFlux = max(ensemble.reactionFluxAllosteric(promisc_rxns_list));   % promiscuous reactions are independent
                end
            else
                if (size(Nelem,2)>1)
                    branchFactor = zeros(size(Nelem,2),1);           
                    for ix = 1:size(Nelem,2)
                        aBranch            = randg(ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:));                    
                        branchFactor(ix,1) = aBranch(1)/sum(aBranch);
                    end
                end
            end
            models(1).rxnParams(activRxnIdx).branchFactor = branchFactor';

            % D. Sample modifier elementary fluxes (positions are given were exp(R)=1)
            modifierElemFlux = models(1).rxnParams(activRxnIdx).modiferElemFlux';
            
            % IV. Calculate rate parameters
            %reactionFlux = ensemble.fluxRef(ensemble.kinActRxns(activRxnIdx));              
            
			
            % VI. Calculate rate parameters
            forwardFlux    = ensemble.forwardFlux{ensemble.kinActRxns(activRxnIdx),strucIdx};
			models(1).rxnParams(activRxnIdx).kineticParams = ...
                calculateKineticParams(reverTemp,forwardFlux,reactionFlux,randomEnzymesR,Nelem,branchFactor,modifierElemFlux);
        end                
    end   

    % Test model consistency
    kineticFxn = str2func(ensemble.kineticFxn{strucIdx});	
    testFlux   = feval(kineticFxn,ones(size(ensemble.freeVars,1),1),models,ensemble.fixedExch(:,1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);
    %disp(testFlux);
    
    % If the model is consistent continue
    if all(abs(testFlux-ensemble.fluxRef)<1e-6)

        % Simulate fluxes
		tolScore      = [];
		condition     = true;
        simulatedFlux = zeros(numel(ensemble.activeRxns),ensemble.numConditions);
        xopt          = zeros(size(x0,1),ensemble.numConditions);
        for ix = 1:ensemble.numConditions

			% % Define anonymous function (comment lines 245-246, and uncomment the lines below to properly simulate)
			% opt.min_objective = @(x) kineticFxn(x,models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);
			% opt.lower_bounds  = lb(:,ix);
            % opt.upper_bounds  = ub(:,ix);

            % % Solve S*v(k,X) = 0; s.t. A*X <= beq, lb < X <ub, with extra constraints (e.g., pool or ratio constraints). Otherwise solve solve S*v(k,X) = 0; s.t. lb < X <ub, with no extra constraints
            % if ~isempty(ensemble.poolConst)
				% for jx = 1:numel(ensemble.poolConst)					
                    % opt.fc{1,2*jx-1} = (@(x) poolConstraintFxn(x,[A_opt{jx},zeros(1,numel(ensemble.freeVars)-numel(ensemble.metsActive))],b{jx}(2*ix-1)));
					% opt.fc{1,2*jx}   = (@(x) poolConstraintFxn(x,[-A_opt{jx},zeros(1,numel(ensemble.freeVars)-numel(ensemble.metsActive))],-b{jx}(2*ix)));					
                % end                
                % opt.fc_tol  = 1e-6*ones(1,2*numel(ensemble.poolConst));                                
            % end
			
			% % FMINCON call
			% [xopt(:,ix),fmin] = fmincon(kineticFxn,x0(:,ix),[],[],[],[],lb(:,ix),ub(:,ix),[],options,models,ensemble.fixedExch,ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},1);		
			
			% % NLOPT call
			% [xopt(:,ix),fmin] = nlopt_optimize(opt,x0(:,ix));
			
			% Simulate fluxes
			xopt(:,ix)          = x0;			% comment this line if you need to simulate
			fmin                = 0;			% comment this line if you need to simulate
            simulatedFlux(:,ix) = feval(kineticFxn,xopt(:,ix),models,ensemble.fixedExch(:,ix+1),ensemble.Sred,ensemble.kinInactRxns,ensemble.subunits{strucIdx},0);             
                   
            % Check mass balance consistency
            if (fmin<massTol)				
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
            tolScore       = max(tolScore);
            acceptanceRate = 1/counter;
            break;
            
        % Erase generated structure if it is not mass balanced
        else
            models = [];
        end
	else
		disp(['There are consistency problems. ModelIdx: ',num2str(strucIdx),', RxnIdxs: ',num2str(find(abs(testFlux-ensemble.fluxRef)>1e-6)),...
		', Error: ',num2str(max(abs(testFlux-ensemble.fluxRef)))]);
	end
end

% Save results and write progress to file
try
	load progress.txt
	progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
	save progress.txt -ascii progress
	save(['tempRejection/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux');
	
% If another worker is writing on the file, wait a small random time
catch
	pause(randi(2)*rand(1));
	load progress.txt
	progress = [progress(1)+1;progress(4)/(progress(1)+1);progress(5)/(progress(1)+1);progress(4)+acceptanceRate;progress(5)+tolScore];
	save progress.txt -ascii progress
	save(['tempRejection/particle_',num2str(progress(1)),'.mat'],'models','strucIdx','xopt','tolScore','simulatedFlux');
end