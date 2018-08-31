function log_particleWeight = computeParticleWeightsVectorizedMGK(models,strucIdx,ensemble,popIdx,prevParams)
%--------------------------------------------------------------------------
% Compute (unnormalized) particle weight
%
% Inputs:       models (structure)
%
% Outputs:      particle unnormalized weight (double)
%--------------------- Pedro Saa 2016 -------------------------------------
% Initialize prior probability (model probability first)
logp_prior = log(1/ensemble.numStruct);															% model structure prior
paramArray = [];																				% initialize parameter array
        
%% 1. Compute prior probability of the proposed particle
%  I. Extract pool constraint information
if ~isempty(ensemble.poolConst)
    for ix = 1:numel(ensemble.poolConst)
        alphaPrior          = ensemble.populations(1).probParams(strucIdx).alphaPoolFactor{ix};
        [unconsPool,U_full] = map2uncons(models(1).poolFactor{ix});	
        logp_prior          = logp_prior + logp_ilrtMKG(unconsPool,alphaPrior,U_full);
        paramArray          = [paramArray,unconsPool];
    end
end                

%  II. Extract gibbs free energy
if ~isempty(ensemble.thermoActive)    
	muPrior    = ensemble.populations(1).probParams(strucIdx).muGibbsFactor;
	sigmaPrior = ensemble.populations(1).probParams(strucIdx).sigmaGibbsFactor;
	gibbsTemp  = log(models(1).gibbsFactor./(1-models(1).gibbsFactor));
	logp_prior = logp_prior + logmvnpdf(gibbsTemp(:)',muPrior(:)',sigmaPrior);
	paramArray = [paramArray,gibbsTemp(:)'];
end

% Loop through all the reactions and initialize hyperparameters
for activRxnIdx = 1:numel(ensemble.kinActRxns)

    % Skip non-enzymatic rxns
	if strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'diffusion')||...
		strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'freeExchange')||...
		strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'fixedExchange')||...
		strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'massAction')
		continue;
	
	% Compute prior probabilities for the rxn-specific parameters
	else		
    
		% II. Compute weights of the enzyme abundances
        [unconsEnz,U_full] = map2uncons(models(1).rxnParams(activRxnIdx).enzymeAbundances);
        alphaPrior         = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances;
		logp_prior         = logp_prior + logp_ilrtMKG(unconsEnz,alphaPrior,U_full);
		paramArray         = [paramArray,unconsEnz];
		
        % III. Compute weights of the reversibilities
        revMatrix 			= ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
        revPrev   			= models(1).rxnParams(activRxnIdx).reversibilities;
		revPrev(revPrev==0) = [];																					% Remove fixed reversibilites
		[unconsRev,U_full]  = map2uncons(revPrev);
        alphaPrior          = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;		
        logp_prior          = logp_prior + logp_ilrtMKG(unconsRev,alphaPrior,U_full);
		paramArray          = [paramArray,unconsRev];
		
        % IV. Compute weights for branch factors (if any) (here we are using the logit-beta transformation)
        Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
        if (size(Nelem,2)>1)
            for ix = 1:size(Nelem,2)
                branchTemp      = models(1).rxnParams(activRxnIdx).branchFactor(ix);
				branchTempTrans = log(branchTemp/(1-branchTemp));
                betaPrior       = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:);
				logp_prior   	= logp_prior + logp_logitBeta(branchTempTrans,betaPrior);
				paramArray      = [paramArray,branchTempTrans];
            end
        end

        % V. Compute weights for mod elem flux (if any) (here we are using the logit-beta transformation)
        if (size(revMatrix,1)==1)&&any(revMatrix==0)
            for ix = 1:sum(revMatrix==0)
                modElemTemp      = models(1).rxnParams(activRxnIdx).modiferElemFlux(ix);
				modElemTempTrans = log(modElemTemp/(1-modElemTemp));
                betaPrior        = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:);
				logp_prior    	 = logp_prior + logp_logitBeta(modElemTempTrans,betaPrior);
				paramArray       = [paramArray,modElemTempTrans];
            end
        end

        % VI. Check whether the current reaction is allosteric. If so, compute prior
        if ensemble.allosteric{strucIdx}(activRxnIdx)			
			allostericTemp      = models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).allostericFactors;
			allostericTempTrans = log(allostericTemp./(1-allostericTemp));
			muPrior             = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).muAllostericParams;
			sigmaPrior          = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).sigmaAllostericParams;
			logp_prior       	= logp_prior + logmvnpdf(allostericTempTrans(:)',muPrior(:)',sigmaPrior);
			paramArray          = [paramArray,allostericTempTrans(:)'];
		end
	end
end

%% 2. Compute importance weight of the proposed particle
% Initialize denominator weight with model probability
logp_weight_struct = log(ensemble.populations(popIdx-1).structureWeights(strucIdx));								% Based on the expression of Toni et al. 2010

% Compute transition probability for the current proposed particle (vectorized version)
sigma       = ensemble.populations(popIdx).probParams(strucIdx).globalKernelParam;									% Extract model-specific cov matrix
logp_weight = log(ensemble.populations(popIdx-1).weights(ensemble.populations(popIdx-1).strucIdx==strucIdx));
logp_transition_kernel = logmvnpdf(paramArray,prevParams(strucIdx).paramArray',sigma);
		
% Add weight contribution of the current active particle
logp_importance_proposal = logp_weight_struct + logsumexp(logp_weight + logp_transition_kernel);

% Compute weight of current inactive particle
log_particleWeight = logp_prior - logp_importance_proposal;