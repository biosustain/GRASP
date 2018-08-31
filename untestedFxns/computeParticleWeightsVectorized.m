function log_particleWeight = computeParticleWeightsVectorized(models,strucIdx,ensemble,popIdx,prevParams)
%--------------------------------------------------------------------------
% Compute (unnormalized) particle weight
%
% Inputs:       models (structure)
%
% Outputs:      particle unnormalized weight (double)
%--------------------- Pedro Saa 2016 -------------------------------------
% Determine alive particles for weight computation
aliveParticles = find(ensemble.populations(popIdx-1).strucIdx==strucIdx);

% Initialize prior probability (model probability first)
log_priorProb = log(1/ensemble.numStruct);
        
%% Compute prior probability of the proposed particle
% 1. Compute weights for pool factors (if any)
if ~isempty(ensemble.poolConst)            
    alphaPrior    = ensemble.populations(1).probParams(strucIdx).alphaPoolFactor;
    log_priorProb = log_priorProb + log(probDirichlet(alphaPrior,models(1).poolFactor,1));
end                

% 2. Compute MLE for the gibbs factors
if ~isempty(ensemble.thermoActive)    
	muPrior       = ensemble.populations(1).probParams(strucIdx).muGibbsFactor;
	sigmaPrior    = ensemble.populations(1).probParams(strucIdx).sigmaGibbsFactor;
	gibbsTemp     = log(models(1).gibbsFactor./(1-models(1).gibbsFactor));
	log_priorProb = log_priorProb + log(mvnpdf(gibbsTemp(:)',muPrior(:)',sigmaPrior));
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
    
		% I. Compute weights of the enzyme abundances
        enzRTemp      = models(1).rxnParams(activRxnIdx).enzymeAbundances;
        alphaPrior    = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaEnzymeAbundances;
        log_priorProb = log_priorProb + log(probDirichlet(alphaPrior,enzRTemp,1));

        % II. Compute weights of the reversibilities
        revMatrix     = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
        revPrev       = models(1).rxnParams(activRxnIdx).reversibilities;
		revPrev(revPrev==0) = [];																					% Remove fixed reversibilites
        alphaPrior    = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;		
        log_priorProb = log_priorProb + log(probDirichlet(alphaPrior,revPrev,1));
		
        % III. Compute weights for branch factors (if any)
        Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
        if (size(Nelem,2)>1)
            for ix = 1:size(Nelem,2)
                branchTemp    = models(1).rxnParams(activRxnIdx).branchFactor(ix);
                betaPrior     = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaBranchFactor(ix,:);                        
                log_priorProb = log_priorProb + log(probDirichlet(betaPrior,[branchTemp,1-branchTemp],1));
            end
        end

        % 6. Compute weights for mod elem flux (if any)
        if (size(revMatrix,1)==1)&&any(revMatrix==0)
            for ix = 1:sum(revMatrix==0)
                modElemTemp = models(1).rxnParams(activRxnIdx).modiferElemFlux(ix);
                betaPrior   = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).betaModiferElemFlux(ix,:);
                log_priorProb = log_priorProb + log(probDirichlet(betaPrior,[modElemTemp,1-modElemTemp],1));
            end
        end

        % Check whether the current reaction is allosteric
        if ensemble.allosteric{strucIdx}(activRxnIdx)			
			allostericTemp = models(1).rxnParams(ensemble.kinActRxns(activRxnIdx)).allostericFactors;
			allostericTemp = log(allostericTemp./(1-allostericTemp));
			muPrior       = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).muAllostericParams;
			sigmaPrior    = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).sigmaAllostericParams;
			log_priorProb = log_priorProb + log(mvnpdf(allostericTemp(:)',muPrior(:)',sigmaPrior));
		end
	end
end

%% Compute importance weight of the proposed particle
% Initialize denominator weight with model probability
log_importance_weight = 0;    
log_importance_weight_struct = log(ensemble.populations(popIdx-1).structureWeights(strucIdx)); % Based on the expression of Toni et al. 2010

% Compute transition probability for the current proposed particle (vectorized version)
% 1. Compute weights for pool factors (if any)
if ~isempty(ensemble.poolConst)
    poolFactorPrev = prevParams(strucIdx).poolFactor;
    KpoolFactor    = ensemble.populations(popIdx).probParams(strucIdx).KpoolFactor;
    poolTemp       = models(1).poolFactor;                
    log_importance_weight = log_importance_weight + log(probDirichlet(KpoolFactor*poolFactorPrev+1,poolTemp,4));
end
	
% 2. Compute weights for the gibbs factors
if ~isempty(ensemble.thermoActive)
	gibbsFactorPrev = prevParams(strucIdx).gibbsFactor;		
	gibbsFactorPrev = log(gibbsFactorPrev./(1-gibbsFactorPrev));
	sigmaGibbs      = ensemble.populations(popIdx).probParams(strucIdx).sigmaGibbsFactor;
	gibbsTemp       = log(models(1).gibbsFactor./(1-models(1).gibbsFactor));
	log_importance_weight = log_importance_weight + log(mvnpdf(gibbsFactorPrev,gibbsTemp(:)',sigmaGibbs));
end
     
% Loop thorugh all the reactions and initialize hyperparameters
for activRxnIdx = 1:numel(ensemble.kinActRxns)
            
	% Skip non-enzymatic rxns
	if strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'diffusion')||...
		strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'freeExchange')||...
		strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'fixedExchange')||...
		strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'massAction')
		continue;
	
	% Compute importance weight for the rxn-specific parameters
	else
					
		% 1. Compute weights of the enzyme abundances
		enzRPrev   = prevParams(strucIdx).enzymeAbundances{activRxnIdx};
		KenzymesR  = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KenzymeAbundances;
		enzRTemp   = models(1).rxnParams(activRxnIdx).enzymeAbundances;
		log_importance_weight = log_importance_weight + log(probDirichlet(KenzymesR*enzRPrev+1,enzRTemp,4));

		% 2. Compute weights of the reversibilities
		revMatrix        = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
		revPrev          = prevParams(strucIdx).reversibilities{activRxnIdx};
		Kreversibilities = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).Kreversibilities;
		revTemp          = models(1).rxnParams(activRxnIdx).reversibilities;
		revTemp(revTemp==0)   = [];																					% Remove invalid reversibilities
		log_importance_weight = log_importance_weight + log(probDirichlet(Kreversibilities*revPrev+1,revTemp,4));

		% 3. Compute weights for branch factors (if any)
		Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
		if (size(Nelem,2)>1)
			for ix = 1:size(Nelem,2)
				branchPrev    = prevParams(strucIdx).branchFactor{activRxnIdx}(:,ix);
				KbranchFactor = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KbranchFactor(ix);				
				branchTemp    = models(1).rxnParams(activRxnIdx).branchFactor(ix);
				log_importance_weight = log_importance_weight + log(probDirichlet(KbranchFactor*[branchPrev,1-branchPrev]+1,[branchTemp,1-branchTemp],4));
			end
		end

		% 4. Compute weights for mod elem flux (if any)
		if (size(revMatrix,1)==1)&&any(revMatrix==0)
			for ix = 1:sum(revMatrix==0)
				modElemPrev  = prevParams(strucIdx).modiferElemFlux{activRxnIdx}(:,ix);
				KmodElemTemp = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KmodiferElemFlux(ix);
				elemTemp     = models(1).rxnParams(activRxnIdx).modiferElemFlux(ix);
				log_importance_weight = log_importance_weight + log(probDirichlet(KmodElemTemp*[modElemPrev,1-modElemPrev]+1,[elemTemp,1-elemTemp],4));
			end
		end

		% Check whether the current reaction is allosteric
		if ensemble.allosteric{strucIdx}(activRxnIdx)            
			allostericPrev        = prevParams(strucIdx).allostericFactors{activRxnIdx};
			allostericPrev        = log(allostericPrev./(1-allostericPrev));																			% apply logit transformation
			sigmaAllosteric       = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).sigmaAllostericParams;
			allostericTemp        = log(models(1).rxnParams(activRxnIdx).allostericFactors./(1-models(1).rxnParams(activRxnIdx).allostericFactors));	% apply logit transformation
			log_importance_weight = log_importance_weight + log(mvnpdf(allostericPrev,allostericTemp(:)',sigmaAllosteric));            
		end
	end		
end
			
% Add weight contribution of the current active particle
log_importance_weight = logsumexp(log_importance_weight + log(ensemble.populations(popIdx-1).weights(aliveParticles)));

% Compute weight of current inactive particle
log_particleWeight = log_priorProb - log_importance_weight - log_importance_weight_struct;