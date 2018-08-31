function ensemble = computeKernelParams(ensemble,popIdx,method)
%--------------------------------------------------------------------------
% Compute kernel parameters 
%
% Inputs:       ensemble (structure) , popIdx (double)
%
% Outputs:      ensemble kernel parameters updated (structure)
%--------------------- Pedro Saa 2016 -------------------------------------
% Loop through all the structures and rxns in the model
for strucIdx = 1:ensemble.numStruct

    % Extract alive particles from the current model structure
    idxParticleStructure = find((ensemble.populations(popIdx-1).strucIdx==strucIdx)&(ensemble.populations(popIdx-1).weights~=0));
    weights = log(ensemble.populations(popIdx-1).weights(idxParticleStructure)) - max(log(ensemble.populations(popIdx-1).weights(idxParticleStructure)));		% Normalize weigths
	weights = exp(weights)/sum(exp(weights));
	
    % Skip structure with no valid particles
    if isempty(idxParticleStructure); continue; end;        
    
    % 1. Compute MLE for the pool factors (if any)
    if ~isempty(ensemble.poolConst)
        
        % Retrieve information from the previous population
        poolFactorPrev = zeros(numel(idxParticleStructure),sum(ensemble.poolConst(1,1:end-1)~=0));
        for ix = 1:numel(idxParticleStructure)
            poolFactorPrev(ix,:) = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).poolFactor;
        end
        
        % Compute kernel parameter based on previous population
        ensemble.populations(popIdx).probParams(strucIdx).KpoolFactor = getKernelParam(poolFactorPrev,method,weights);
    end
    
	% 2. Compute MLE for the gibbs factors
	if ~isempty(ensemble.thermoActive)
        gibbsFactorPrev = zeros(numel(idxParticleStructure),numel(ensemble.thermoActive));
        for ix = 1:numel(idxParticleStructure)
            gibbsFactorPrev(ix,:) = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).gibbsFactor;
        end
                
		% Compute kernel parameter based on previous population
		gibbsFactorPrev = log(gibbsFactorPrev./(1-gibbsFactorPrev));			% apply logit transformation
		muX   = weights'*gibbsFactorPrev;
		diffX = gibbsFactorPrev-muX(ones(size(gibbsFactorPrev,1),1),:);
		diffW = weights(:,ones(1,size(gibbsFactorPrev,2))).*diffX;
		covX  = diffW'*diffX/(1-sum(weights.^2));		
		ensemble.populations(popIdx).probParams(strucIdx).muGibbsFactor    = muX;
		ensemble.populations(popIdx).probParams(strucIdx).sigmaGibbsFactor = 2*covX;	% twice the weighted sample covariance
    end	
	
    % II. Loop thorugh all the reactions and initialize hyperparameters
    for activRxnIdx = 1:numel(ensemble.kinActRxns)
        
        % Compute kernel parameters for all the rxns except the non-enzymatic
		if strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'diffusion')||...
			strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'freeExchange')||...
			strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'fixedExchange')||...
			strcmp(ensemble.rxnMechanisms{strucIdx}{ensemble.kinActRxns(activRxnIdx)},'massAction')
			continue;

		% Compute probability parameters for the enzymatic mechanisms
        else
                                                            
            % 1. Retrieve EnzR information from the previous population            
            enzRPrev = zeros(numel(idxParticleStructure),max(ensemble.forwardFlux{ensemble.kinActRxns(activRxnIdx),strucIdx}(:)));
            for ix = 1:numel(idxParticleStructure)
                enzRPrev(ix,:) = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).rxnParams(ensemble.kinActRxns(activRxnIdx)).enzymeAbundances;
            end
            
            % Compute kernel parameter based on previous population
            ensemble.populations(popIdx).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).KenzymeAbundances = getKernelParam(enzRPrev,method,weights);
                        
            % 2. Retrieve Rev information from the previous population
            revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
            revPrev   = zeros(numel(idxParticleStructure),max([all(size(revMatrix,1)==1)*sum(revMatrix~=0),(size(revMatrix,1)>1)*size(revMatrix,2)]));
            for ix = 1:numel(idxParticleStructure)
                revTemp       = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).rxnParams(ensemble.kinActRxns(activRxnIdx)).reversibilities;
                revPrev(ix,:) = revTemp(revTemp~=0);
            end
            
            % Compute kernel parameter based on previous population
            ensemble.populations(popIdx).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).Kreversibilities = getKernelParam(revPrev,method,weights);
            
            % 3. Retrieve branch information from the previous population (if any)
            Nelem = ensemble.Nelem{ensemble.kinActRxns(ensemble.kinActRxns(activRxnIdx)),strucIdx};
            
            % If there is branch information available
            if (size(Nelem,2)>1)
                branchPrev = zeros(numel(idxParticleStructure),size(Nelem,2));
                for ix = 1:numel(idxParticleStructure)
                    branchPrev(ix,:) = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).rxnParams(ensemble.kinActRxns(activRxnIdx)).branchFactor;
                end
                for ix = 1:size(branchPrev,2)
                    
                    % Compute kernel parameter based on previous population
                    ensemble.populations(popIdx).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).KbranchFactor(ix,1) = getKernelParam([branchPrev(:,ix),1-branchPrev(:,ix)],method,weights);
                end
            end
            
            % 4. Retrieve information about modifiers flux (if any)
            if (size(revMatrix,1)==1)&&any(revMatrix==0)
                modElemFluxPrev = zeros(numel(idxParticleStructure),sum(revMatrix==0));
                for ix = 1:numel(idxParticleStructure)
                    modElemFluxPrev(ix,:) = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).rxnParams(ensemble.kinActRxns(activRxnIdx)).modiferElemFlux;
                end
                
                % Solve Beta parameters for the branch factors
                for ix = 1:sum(revMatrix==0)
                    
                    % Compute kernel parameter based on previous population
                    ensemble.populations(popIdx).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).KmodiferElemFlux(ix,1) = getKernelParam([modElemFluxPrev(:,ix),1-modElemFluxPrev(:,ix)],method,weights);
                end
            end
            
            % Check whether the current reaction is allosteric
            if ensemble.allosteric{strucIdx}(ensemble.kinActRxns(activRxnIdx))               
				
				% Compute probability parameters for the allosteric mxn				
				allostericPrev = zeros(numel(idxParticleStructure),numel(ensemble.populations(1).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).muAllostericParams));
                for ix = 1:numel(idxParticleStructure)
					allostericPrev(ix,:) = ensemble.populations(popIdx-1).models(idxParticleStructure(ix)).rxnParams(ensemble.kinActRxns(activRxnIdx)).allostericFactors;
                end               
				
				% Compute kernel parameter based on previous population
				allostericPrev = log(allostericPrev./(1-allostericPrev));			% apply logit transformation
				muX   = weights'*allostericPrev;
				diffX = allostericPrev-muX(ones(size(allostericPrev,1),1),:);
				diffW = weights(:,ones(1,size(allostericPrev,2))).*diffX;
				covX  = diffW'*diffX/(1-sum(weights.^2));
				ensemble.populations(popIdx).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).muAllostericParams	   = muX;		
				ensemble.populations(popIdx).probParams(strucIdx).rxnParams(ensemble.kinActRxns(activRxnIdx)).sigmaAllostericParams = 2*covX; 	% twice the weighted sample covariance
            end
        end
    end
end
disp('Kernel parameters updated.');