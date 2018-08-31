function log_particleWeight = computeParticleWeights(models,strucIdx,ensemble,popIdx)
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
        revPrevValid  = revPrev(revPrev~=0);                                                                     % remove fixed reversibilities
        alphaPrior    = ensemble.populations(1).probParams(strucIdx).rxnParams(activRxnIdx).alphaReversibilities;
        log_priorProb = log_priorProb + log(probDirichlet(alphaPrior,revPrevValid,1));

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
log_importance_weight = 0;    % Based on the expression of Toni et al. 2010, otherwise = log(ensemble.populations(popIdx-1).structureWeights(strucIdx))
        
% Compute transition probability for the current proposed particle
for jx = 1:numel(aliveParticles)
            
	% 1. Compute weights for pool factors (if any)
    if ~isempty(ensemble.poolConst)
        poolFactorPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).poolFactor;
        KpoolFactor    = ensemble.populations(popIdx).probParams(strucIdx).KpoolFactor;
        poolTemp       = models(1).poolFactor;                
        log_importance_weight = log_importance_weight + log(probDirichlet(KpoolFactor*poolFactorPrev+1,poolTemp,1));
    end
	
	% 2. Compute weights for the gibbs factors
	if ~isempty(ensemble.thermoActive)
		gibbsFactorPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).gibbsFactor;		
		gibbsFactorPrev = log(gibbsFactorPrev./(1-gibbsFactorPrev));
		sigmaGibbs      = ensemble.populations(popIdx).probParams(strucIdx).sigmaGibbsFactor;
		gibbsTemp       = log(models(1).gibbsFactor./(1-models(1).gibbsFactor));
		log_importance_weight = log_importance_weight + log(mvnpdf(gibbsFactorPrev(:)',gibbsTemp(:)',sigmaGibbs));
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
			enzRPrev   = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).enzymeAbundances;
			KenzymesR  = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KenzymeAbundances;
			enzRTemp   = models(1).rxnParams(activRxnIdx).enzymeAbundances;
			log_importance_weight = log_importance_weight + log(probDirichlet(KenzymesR*enzRPrev+1,enzRTemp,1));

			% 2. Compute weights of the reversibilities
			revMatrix        = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
			revPrev          = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).reversibilities;
			Kreversibilities = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).Kreversibilities;
			revTemp          = models(1).rxnParams(activRxnIdx).reversibilities;
			log_importance_weight = log_importance_weight + log(probDirichlet(Kreversibilities*revPrev+1,revTemp,1));

			% 3. Compute weights for branch factors (if any)
			Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
			if (size(Nelem,2)>1)
				for ix = 1:size(Nelem,2)
					branchPrev    = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).branchFactor(ix);
					KbranchFactor = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KbranchFactor(ix);
					branchTemp    = models(1).rxnParams(activRxnIdx).branchFactor(ix);
					log_importance_weight = log_importance_weight + log(probDirichlet(KbranchFactor*[branchPrev,1-branchPrev]+1,[branchTemp,1-branchTemp],1));
				end
			end

			% 4. Compute weights for mod elem flux (if any)
			if (size(revMatrix,1)==1)&&any(revMatrix==0)
				for ix = 1:sum(revMatrix==0)
					modElemPrev  = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).modiferElemFlux(ix);
					KmodElemTemp = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).KmodiferElemFlux(ix);
					elemTemp     = models(1).rxnParams(activRxnIdx).modiferElemFlux(ix);
					log_importance_weight = log_importance_weight + log(probDirichlet(KmodElemTemp*[modElemPrev,1-modElemPrev]+1,[elemTemp,1-elemTemp],1));
				end
			end

			% Check whether the current reaction is allosteric
			if ensemble.allosteric{strucIdx}(activRxnIdx)            
				allostericPrev        = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).allostericFactors;			
				allostericPrev        = log(allostericPrev./(1-allostericPrev));																			% apply logit transformation
				sigmaAllosteric       = ensemble.populations(popIdx).probParams(strucIdx).rxnParams(activRxnIdx).sigmaAllostericParams;
				allostericTemp        = log(models(1).rxnParams(activRxnIdx).allostericFactors./(1-models(1).rxnParams(activRxnIdx).allostericFactors));	% apply logit transformation
				log_importance_weight = log_importance_weight + log(mvnpdf(allostericPrev(:)',allostericTemp(:)',sigmaAllosteric));            
			end
		end
	end
			
	% Add weight contribution of the current active particle
	log_importance_weight = log_importance_weight + log(ensemble.populations(popIdx-1).weights(aliveParticles(jx)));
end       

% Compute weight of current inactive particle
log_particleWeight = log_priorProb-log_importance_weight;