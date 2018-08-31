function prevParams = getPopulationParams(ensemble,popIdx)
%--------------------------------------------------------------------------
% Assemble parameters from the previous population
%
% Inputs:       ensemble (structure), popIdx (double)
%
% Outputs:      preParams (structure array)
%--------------------- Pedro Saa 2017 -------------------------------------
% Initialize structure array
for strucIdx = 1:ensemble.numStruct
	prevParams(strucIdx).poolFactor  = [];
	prevParams(strucIdx).gibbsFactor = [];
	prevParams(strucIdx).enzymeAbundances{numel(ensemble.kinActRxns),1}  = [];
	prevParams(strucIdx).reversibilities{numel(ensemble.kinActRxns),1}   = [];
	prevParams(strucIdx).branchFactor{numel(ensemble.kinActRxns),1}      = [];
	prevParams(strucIdx).modiferElemFlux{numel(ensemble.kinActRxns),1}   = [];
	prevParams(strucIdx).allostericFactors{numel(ensemble.kinActRxns),1} = [];
end

% Loop through all the structures in the model
for strucIdx = 1:ensemble.numStruct

	% Check that the structure is valid
	if (ensemble.populations(popIdx-1).structureWeights(strucIdx)==0)
		continue;

	% Proceed for this model structure
	else	
		aliveParticles = find(ensemble.populations(popIdx-1).strucIdx==strucIdx);
		
		% Loop through all the alive particles
		for jx = 1:numel(aliveParticles)
            
			% 1. Compute weights for pool factors (if any)
			if ~isempty(ensemble.poolConst)
				poolFactorPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).poolFactor;
				prevParams(strucIdx).poolFactor = [prevParams(strucIdx).poolFactor;poolFactorPrev(:)'];
			end
	
			% 2. Compute weights for the gibbs factors
			if ~isempty(ensemble.thermoActive)
				gibbsFactorPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).gibbsFactor;
				prevParams(strucIdx).gibbsFactor = [prevParams(strucIdx).gibbsFactor;gibbsFactorPrev(:)'];
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
					enzRPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).enzymeAbundances;
					prevParams(strucIdx).enzymeAbundances{activRxnIdx} = [prevParams(strucIdx).enzymeAbundances{activRxnIdx};enzRPrev(:)'];
				
					% 2. Compute weights of the reversibilities
					revMatrix = ensemble.revMatrix{ensemble.kinActRxns(activRxnIdx),strucIdx};
					revPrev   = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).reversibilities;
					revPrev(revPrev==0) = [];																						% Remove fixed reversibilities
					prevParams(strucIdx).reversibilities{activRxnIdx} = [prevParams(strucIdx).reversibilities{activRxnIdx};revPrev(:)'];			
					
					% 3. Compute weights for branch factors (if any)
					Nelem = ensemble.Nelem{ensemble.kinActRxns(activRxnIdx),strucIdx};
					if (size(Nelem,2)>1)						
						for ix = 1:size(Nelem,2)
							branchPrev(1,ix) = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).branchFactor(ix);					
						end
						prevParams(strucIdx).branchFactor{activRxnIdx} = [prevParams(strucIdx).branchFactor{activRxnIdx};branchPrev];
					end
																				
					% 4. Compute weights for mod elem flux (if any)
					if (size(revMatrix,1)==1)&&any(revMatrix==0)
						for ix = 1:sum(revMatrix==0)
							modElemPrev(1,ix) = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).modiferElemFlux(ix);
						end
						prevParams(strucIdx).modiferElemFlux{activRxnIdx} = [prevParams(strucIdx).modiferElemFlux{activRxnIdx};modElemPrev];
					end

					% Check whether the current reaction is allosteric
					if ensemble.allosteric{strucIdx}(activRxnIdx)            
						allostericPrev = ensemble.populations(popIdx-1).models(aliveParticles(jx)).rxnParams(activRxnIdx).allostericFactors;
						prevParams(strucIdx).allostericFactors{activRxnIdx} = [prevParams(strucIdx).allostericFactors{activRxnIdx};allostericPrev(:)'];
					end
				end
			end
		end			
	end       
end